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
#let CMO = math.op("CMO")
#let PLD = math.op("PLD")
#let BR = math.op("BR")
#let Rev = math.op("Rev")
#let softmax = math.op("softmax")
#let pos = math.op("pos")
#let ind = math.op("1")

#let keybox(title, body) = block(
  width: 100%,
  breakable: false,
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
  width: 100%,
  breakable: false,
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
  width: 100%,
  breakable: false,
  inset: 8pt,
  stroke: 0.68pt + rgb("2f5f7f"),
  radius: 5pt,
  fill: rgb("f4f9fc"),
)[
  #strong(title)
  #v(0.26em)
  #body
]

#let fig-note(body) = block(
  width: 100%,
  breakable: false,
  inset: 7pt,
  stroke: 0.48pt + rgb("8796a3"),
  radius: 4pt,
  fill: rgb("f7f9fb"),
)[
  #text(size: 8.65pt, fill: rgb("34434f"))[#body]
]


#let imgfig(title, path, note: none, width: 100%) = figure(
  block(width: 100%)[
    #align(center)[#image(path, width: width)]
    #if note != none [
      #v(0.20em)
      #fig-note(note)
    ]
  ],
  caption: title,
)

#let nbref(body) = block(
  width: 100%,
  inset: (left: 7pt, right: 7pt, top: 4pt, bottom: 4pt),
  stroke: (left: 2pt + rgb("2f5f7f")),
  fill: rgb("f3f8fb"),
  radius: 2pt,
)[
  #text(size: 8.25pt, fill: rgb("244b64"))[
    *Onde está no notebook:* #body
  ]
]

#align(center)[
  #text(17.2pt, weight: "bold")[MFG - Renewable Energy]

  #v(0.24em)
  #text(13.0pt, weight: "semibold", fill: rgb("2f5f7f"))[
    Relatório final
  ]

  #v(0.52em)
  #text(10.1pt)[João Felipe Vilas Boas]

  #v(0.18em)
  #text(8.8pt, fill: gray.darken(20%))[
    Modelo e resultados dos notebooks 07 e 08
  ]
]

#v(0.80em)

#block(inset: 9pt, stroke: 0.65pt + gray, radius: 5pt, fill: luma(250))[
  *Resumo.* Este relatório desenvolve uma estrutura empírico-matemática para
  estudar como a expansão de fontes renováveis intermitentes altera a carga
  residual, o despacho hidrotérmico, o preço-sombra locacional, o valor da água e
  os incentivos de investimento no Sistema Interligado Nacional. O painel contém 35.040 observações
  horárias de 2025, quatro subsistemas e CMO sem lacunas após a junção documentada.
  O clearing comum possui curtailment solar e eólico explícitos, térmica em blocos,
  déficit, rede, reservatório em energia equivalente e fechamento cíclico em cada
  dia representativo. Sobre esse operador são construídos um laboratório finito de
  dois produtores, um benchmark estático de entrada, um planejador social e um
  MFG dinâmico locacional. O benchmark estático e o MFG dinâmico satisfazem as
  tolerâncias com resíduos calculados antes da relaxação. A trajetória dinâmica
  conserva massa, associa cada política a uma equação forward e mostra a evolução
  conjunta da distribuição, da localização e da capacidade. Os números são
  condicionais às proxies de custo, hidrologia, rede e capacidade; a contribuição
  principal é a cadeia coerente entre mecanismo empírico, clearing físico, decisão
  individual e consistência de campo médio.
]

#v(0.55em)
#text(11pt, weight: "semibold", fill: rgb("2f5f7f"))[
  Guia de rastreabilidade: relatório $arrow.l.r$ notebooks
]

As referências abaixo usam primeiro a seção e o nome do bloco, que permanecem
estáveis. O número entre colchetes é o contador de execução da versão entregue e
pode mudar quando o notebook for reexecutado. Assim, por exemplo, “07 §6.2,
`solve_mfg`, `[17]`” significa: notebook 07, Seção 6.2, célula que define
`solve_mfg`, atualmente executada como `[17]`.

#text(size: 7.7pt)[
#table(
  columns: (1.18fr, 2.25fr, 1.68fr),
  inset: 3.5pt,
  align: left,
  [Objeto no relatório], [Implementação no notebook 07], [Figura/auditoria no notebook 08],
  [Painel e CMO], [§1 e §1.1, células `[3]`--`[6]`; junção por hora/subsistema e cache consolidado.], [§1, célula `[2]`, preço × carga residual.],
  [Capacidades e custos], [§2, célula `[7]`; `q99`, `K_SOLAR_OBS` e `q_type`.], [§9, célula `[10]`, comparação privado × social.],
  [Dias representativos], [§3, célula `[8]`; `build_representative_days`.], [§11, célula `[12]`, pesos e homogeneidade.],
  [Rede e clearing], [§4.1, `[10]`, `EDGES/FMAX`; §4.2, `[11]`, `DispatchEngineCVX`.], [§10, célula `[11]`, topologia.],
  [Reservatório e KKT], [§4.2--4.3, `[11]`--`[12]`; dinâmica, duais e `reservoir_audit.csv`.], [§5, célula `[6]`, trajetória do estoque.],
  [Hidrologia e vertimento], [§4.5--4.5.1, célula `[13]`; cenários, termos quadráticos e diagnóstico regional.], [§6 e §15, células `[7]` e `[16]`.],
  [$Lambda'$ e canibalização], [§4.6, célula `[14]`; diferença finita de carga e varredura de capacidade.], [§15, célula `[16]`.],
  [Laboratório finito], [§5, célula `[15]`; `SIT_LOC`, incumbente, entrante e `finite_two_player_sensitivity.csv`.], [§7, célula `[8]`.],
  [Benchmark estático], [§6.1--6.4, `[16]`--`[19]`; `best_response`, `solve_mfg` e calibração.], [§8, célula `[9]`, resíduos e capacidade.],
  [Planejador social], [§8, célula `[22]`; `solve_social_planner`.], [§9, célula `[10]`.],
  [MFG dinâmico], [§10, célula `[23]`; `solve_dynamic_bellman`, `dynamic_forward` e `solve_dynamic_mfg`.], [§12--14, células `[13]`--`[15]`.],
  [Tabelas e figuras finais], [§10.1, célula `[24]`; grava os CSVs dinâmicos.], [§16, mapa dos nomes dos arquivos.],
)
]



= Introdução

A expansão solar e eólica muda a pergunta econômica central da operação elétrica.
Não basta comparar energia renovável acumulada com consumo acumulado. É preciso
saber *quando* e *onde* essa energia aparece, se a rede consegue transportá-la, se
o sistema dispõe de potência e rampas nas horas sem sol e como o uso da água hoje
altera a flexibilidade de amanhã.

No Brasil, esses mecanismos aparecem ao mesmo tempo. A solar reduz a carga
residual no meio do dia, a eólica é espacialmente concentrada, a hidro transfere
energia entre horas e os intercâmbios conectam subsistemas com disponibilidades
diferentes. Quando muitos projetos solares produzem no mesmo perfil, a oferta
agregada reduz o valor justamente nas horas em que esses projetos vendem. Essa
realimentação é a canibalização de receita que conecta a parte empírica ao jogo de
campo médio.

A pergunta principal é:

#align(center)[
  #box(inset: 8pt, stroke: 0.65pt + black, radius: 4pt)[
    *Como o aumento da participação de fontes renováveis intermitentes altera o
    despacho hidrotérmico, o preço de curto prazo modelado e os incentivos privados
    de investimento em um mercado locacional de energia?*
  ]
]

A resposta é organizada em cinco camadas causais:

#enum[
  dados horários e fatos estilizados de carga residual, CMO e captura;
][
  clearing locacional com rede, curtailment, déficit e reservatório explícitos;
][
  laboratório finito que isola incumbente e entrante;
][
  investimento estático privado, realocação com total fixo e planejador social;
][
  MFG dinâmico com Bellman, política probabilística, equação forward, clearing e
  ponto fixo sobre a trajetória inteira.
]


= Dados, notação e compressão temporal

== Índices, estados e unidades

O tempo operacional é horário, $t=0,dots,23$. Os dias representativos são
$d in cal(D)$, os subsistemas são

$ cal(L) = ("SE", "S", "NE", "N"), $

os blocos térmicos são $b in cal(B)^"th"$ e os tipos solares são
$theta in ("high", "low")$. No MFG dinâmico, $y=0,dots,H$ denota o estágio de
decisão e não a hora do despacho.

A geração térmica do SIN reúne dezenas de usinas
com custos variáveis muito diferentes — de ciclos combinados eficientes a óleo e
diesel de ponta. Modelar cada máquina exigiria dados de acionamento, rampa mínima e
tempo de partida que não estão disponíveis nesta escala. Em vez disso, agregamos
essa oferta em um pequeno conjunto de *blocos* $b in cal(B)^"th"$, cada um com
custo marginal constante $c_(ell,b)^"th"$ e potência máxima $overline(N)_(ell,b)$.
Um bloco é, portanto, um degrau de uma curva de oferta térmica: um intervalo de
potência que pode ser produzido a um mesmo custo unitário. Ordenados do mais barato
ao mais caro, os blocos formam uma aproximação por partes da curva de mérito: o
clearing só aciona um bloco caro depois de esgotar todos os mais baratos, que é
exatamente a lógica do despacho por ordem de mérito.

Eles são o mecanismo por trás de três resultados
centrais do relatório. Primeiro, o preço-sombra: nas horas em que a térmica é a
tecnologia marginal, o bloco acionado no topo da pilha fixa $lambda$, de modo que a
condição de primeira ordem $lambda_t = c_b^"th" + dots$ da seção do clearing liga
diretamente preço a mérito térmico. Segundo, a adequação: solar não
produzem nas horas sem sol nem sob rampa do entardecer, e são os blocos térmicos
(com a hidro) que preenchem essa carga residual; sem uma pilha térmica com potência
suficiente, o problema recorre ao déficit penalizado. Terceiro, a canibalização: ao
crescer, a solar desloca os blocos térmicos *mais baratos* nas horas de sol,
empurrando a tecnologia marginal — e o preço — para baixo justamente quando os
projetos solares vendem. Sem uma representação explícita do mérito térmico, esse
canal de preço teria de ser imposto por uma curva arbitrária.

Nem todo projeto solar é igual. Em vez de tratar cada empreendimento
individualmente, resumimos a heterogeneidade persistente em dois *tipos*
$theta in ("high", "low")$. Aqui $theta$ não é apenas a irradiação do sítio — esse
efeito locacional já entra no perfil $a_(d,ell,t)^s$ —, mas uma qualidade composta do
projeto e do desenvolvedor: eficiência e disponibilidade dos equipamentos, perdas,
qualidade de execução, escala de compra e acesso a financiamento. Essa leitura torna
factível que as vantagens sejam positivamente correlacionadas, sem impor um
trade-off mecânico: projetos de melhor qualidade podem simultaneamente entregar mais
energia útil por MW e enfrentar menor custo anualizado de capital e operação.

No baseline, o tipo `high` tem produtividade relativa $xi_"high"=1,08$ e
multiplicador de custo 0,94; o tipo `low` tem $xi_"low"=0,94$ e multiplicador 1,08.
As massas iniciais são 45% `high` e 55% `low`. Portanto, `high` é deliberadamente
um tipo dominante nas duas dimensões, e não um ponto estimado de uma fronteira de
custo-desempenho. Essa ordenação é uma hipótese estilizada para representar
qualidade persistente e testar como heterogeneidade afeta seleção e localização; os
valores e as massas não foram estimados com microdados de projetos. Uma formulação
com trade-off exigiria tipos cruzados — por exemplo, alta produtividade com custo
alto — e constitui uma sensibilidade ainda pendente.

No benchmark estático, os tipos separam a decisão de
investimento por qualidade de recurso: a entrada livre e o custo anualizado
$q_(ell,theta)$ dependem de $theta$, e a distribuição por tipo determina quanto da
capacidade é rentável em cada localização. No MFG dinâmico, quem permanece fixo é
$theta$, o tipo solar `high` ou `low`; $b$ não é um tipo, mas o bin de capacidade e
pode mudar a cada estágio. Assim, no estado $x=(ell,theta,b)$, o agente pode migrar
em $ell$ e aumentar ou reduzir $b$, enquanto conserva $theta$. Essa permanência é
uma hipótese de modelagem: interpreta $theta$ como uma característica persistente do
projeto ou investidor, e não como uma qualidade do sítio que seria sorteada novamente
ao migrar. Reduzir tudo a um tipo médio apagaria essa heterogeneidade e daria uma
resposta de investimento artificialmente uniforme.

#keybox([Unidades e objetos que exigem cuidado])[
  Capacidade e potência são medidas em MW; energia horária em MWh; estoques e
  afluências em MWh-equivalentes; receitas anuais por unidade de capacidade em
  R\$/MW/ano; e $lambda$ e CMO em R\$/MWh. A capacidade `p99` é uma proxy
  operacional baseada no percentil 99 da geração observada, não capacidade nominal
  cadastrada. Os pesos dos dias representativos são reescalados para 365 dias e,
  portanto, anualizam o modelo.
]

== Fontes de dados e transformações

Todos os dados observados usados no painel têm referência institucional, arquivo
local e transformação identificáveis. A tabela separa o que entra efetivamente no
baseline do que foi apenas baixado para auditoria ou extensão futura.

#text(size: 7.8pt)[
#table(
  columns: (0.72fr, 1.42fr, 2.10fr, 1.30fr),
  inset: 3.4pt,
  align: left,
  [Objeto], [Fonte original e frequência], [Transformação aplicada], [Uso no relatório],
  [Carga], [ONS, `curva-carga-ho`; `CURVA_CARGA_2025.csv`; horária],
  [Instante convertido para hora, duplicatas agregadas pela média e subsistemas padronizados.],
  [Carga $L$ do painel e do clearing.],
  [Geração por fonte], [ONS, `geracao_usina_2_ho`; 12 arquivos mensais; horária por usina],
  [Valores somados por hora, subsistema e fonte. A solar mantém todas as modalidades do arquivo, inclusive `Pequenas Usinas (MMGD)`.],
  [Perfis observados, carga residual e proxies `p99`.],
  [Afluência equivalente], [Derivada da geração hídrica horária do conjunto anterior.],
  [Soma das 24 horas do medoide, limitada a $24K_ell^"h,p99"$, dividida uniformemente pelas 24 horas.],
  [Entrada $A_(d,ell,t)$ do reservatório equivalente.],
  [CMO], [ONS, `cmo_tm`; `CMO_SEMIHORARIO_2025.csv`; semihorária],
  [Média das duas observações de cada hora; junção por hora e subsistema; preenchimento temporal apenas da lacuna documentada.],
  [Diagnósticos empíricos, fatores de captura e comparação com $lambda$.],
  [Intercâmbio], [ONS, `intercambio_nacional_ho`; `INTERCAMBIO_NACIONAL_2025.csv`; horária],
  [Importação positiva no destino e negativa na origem, agregada como saldo líquido.],
  [Auditoria do painel; o clearing usa a rede reduzida e seus limites, não reproduz o fluxo observado.],
  [Hidrologia física], [ONS, `dados_hidrologicos_di`, ENA e EAR; diária por reservatório.],
  [Níveis e vazões foram preservados para validação e extensões; não alimentam diretamente a afluência do baseline.],
  [Referência física e limitação declarada.],
  [Capacidade cadastrada], [ONS `capacidade-geracao` e ANEEL/SIGA; fotografia cadastral.],
  [Conciliação ONS--SIGA disponível, mas não substitui a proxy operacional na execução corrente.],
  [Sensibilidade futura; não entra no baseline.],
)
]


O recorte solar merece destaque. Como a agregação soma todas as linhas
fotovoltaicas do ONS, $G^s$ e sua capacidade `p99` combinam geração centralizada,
`Pequenas Usinas (Tipo III)`, `TIPO II-B` e MMGD; no arquivo de 2025, a modalidade
MMGD responde por aproximadamente 64,5% da energia solar somada. A carga usada é a
curva $D$ do ONS, sem reconstrução explícita de uma carga bruta adicionando MMGD.
Logo, o baseline não separa mercado centralizado de geração atrás do medidor nem
resolve uma eventual diferença conceitual entre carga medida e carga bruta. Os
resultados solares devem ser lidos como condicionais a essa agregação; separar
centralizada e MMGD e reconciliá-las com a definição da carga é uma validação
empírica necessária antes de uso regulatório.

== Painel horário de 2025

O notebook 07 lê o painel construído, padroniza as colunas e faz
uma junção à esquerda do CMO por `(instante, subsistema)`. O cache de preço não
contém 16 de maio de 2025: são 24 horas ausentes em cada um dos quatro subsistemas,
96 células ou 0,274% do painel. Essas células são preenchidas dentro de cada
subsistema por propagação temporal (`forward fill`, com `backward fill` apenas como
salvaguarda); não há propagação entre regiões. A escolha preserva uma grade completa
para os diagnósticos e para a seleção dos medoides sem inventar diferenças
regionais, mas equivale a supor que o último CMO observado permanece no dia ausente.

#nbref([
  Notebook 07, §1.1, célula `[4]`: procure a função de carregamento do CMO, a
  atribuição `out["CMO_obs"]`, o `merge` por instante/subsistema e
  `groupby("id_subsistema")["CMO_obs"].transform(lambda s: s.ffill().bfill())`.
  A auditoria de cobertura e correlação está em §1.2, célula `[5]`; a gravação do
  painel consolidado está em §1.3, célula `[6]`.
])
É defensável pela pequena participação da lacuna e pela finalidade de continuidade
computacional, não como imputação causal; análises de cauda de preço devem repetir o
cálculo excluindo esse dia. O painel consolidado possui 8.760 horas por
subsistema, de 1º de janeiro a 31 de dezembro de 2025, totalizando 35.040 linhas.

#text(size: 8.8pt)[
#table(
  columns: (0.72fr, 0.88fr, 0.94fr, 0.92fr, 0.92fr, 0.92fr),
  inset: 4pt,
  align: right,
  [Subs.], [Carga média MW], [Carga máxima MW], [Solar média MW], [Eólica média MW], [CMO médio R\$/MWh],
  [SE], [44.225,9], [62.149,9], [5.257,4], [6,2], [216,19],
  [S], [13.804,4], [22.737,4], [1.301,0], [698,8], [226,58],
  [NE], [13.266,9], [16.446,6], [3.341,2], [12.266,5], [150,61],
  [N], [8.311,7], [10.239,2], [583,5], [212,4], [162,79],
)
]

As capacidades efetivas proxy, usadas como referência de calibração e não como
cadastro físico, são calculadas separadamente por fonte e subsistema como
$K_(ell,r)^"p99"=Q_(0,99)({G_(ell,t)^r}_t)$: o percentil 99 das 8.760 somas horárias
de geração. O percentil alto aproxima uma potência efetivamente alcançada sem deixar
que uma única hora extrema determine a escala. Ele não corrige indisponibilidade,
curtailment, expansão ocorrida dentro de 2025 nem a inclusão de MMGD e, portanto,
não deve ser chamado de capacidade de placa. Os valores são:

#text(size: 8.8pt)[
#table(
  columns: (0.75fr, 1fr, 1fr, 1fr, 1fr),
  inset: 4pt,
  align: right,
  [Subs.], [Solar p99 GW], [Eólica p99 GW], [Hidro p99 GW], [Térmica p99 GW],
  [SE], [19,077], [0,016], [44,997], [7,951],
  [S], [5,919], [1,835], [15,243], [2,051],
  [NE], [11,789], [21,842], [6,541], [2,970],
  [N], [2,284], [0,411], [18,715], [3,441],
  [*Total solar*], [*39,069*], [-], [-], [-],
)
]

== Dias representativos

Cada camada de decisão do relatório resolve o
clearing muitas vezes — o MFG dinâmico, por exemplo, o resolve a cada estágio e a
cada iteração do ponto fixo. Rodar os 365 despachos de 24 horas do ano em cada uma
dessas chamadas seria computacionalmente caro. Como muitos dias têm perfis horários
de carga e de geração renovável próximos, adotamos a aproximação de que um conjunto
pequeno preserva os mecanismos qualitativos relevantes — hipótese depois auditada,
mas não garantida. A ideia é
substituir o calendário completo por um punhado de *dias representativos*, cada um
carregando um peso igual ao número de dias reais que ele resume.

Descrevemos cada dia do ano por um vetor de características operacionais agregadas
por subsistema — carga média e de pico, energia solar, eólica e hídrica e CMO médio —
além de marcadores de estação e de calendário, e agrupamos os 365 vetores por
similaridade.
Cada grupo — um *cluster* — reúne dias com comportamento operacional parecido: por
exemplo, dias úteis de inverno com pouca solar, ou fins de semana de verão com forte
geração diurna. Em vez de simular todos os dias do cluster, simulamos apenas um dia
que o represente e anexamos a ele o peso do grupo. O número de clusters é escolhido
entre 4 e 12 combinando três critérios: a *silhueta* (uma medida de quão bem
separados e coesos ficam os grupos), a cobertura das quatro *estações* e a distinção
entre *regimes de calendário* (útil, fim de semana, feriado). A tabela abaixo usa
quatro blocos.

O representante de cada cluster é o seu
*medoide operacional*: o dia real mais próximo do centroide do $k$-médias no espaço
padronizado das características usadas no agrupamento. Ele contrasta com o próprio
*centroide*, que é a média do grupo — um perfil
*sintético* que pode não corresponder a nenhum dia observado e que suaviza
artificialmente rampas e picos ao promediar dias ligeiramente defasados. Como o
clearing depende criticamente da forma horária (o vale solar do meio-dia, a rampa do
entardecer, o horário do pico de carga), preferimos o medoide por três razões: ele
preserva um perfil horário genuinamente observado, sem inventar transições que nunca
ocorreram; é robusto a dias atípicos, que não distorcem um representante real como
distorceriam uma média; e mantém a interpretação física — cada bloco tem uma data
real associada, listada na tabela.

Além da estação, cada dia
recebe um rótulo de *regime de calendário* construído a partir da própria data. Um
dia é *útil* quando é de segunda a sexta e não é feriado: são os dias em que a
atividade industrial e comercial está plena, com carga mais alta e com os dois picos
característicos, o da manhã e o do fim de tarde. Um dia de *fim de semana* é sábado
ou domingo: a carga cai e o perfil se achata, porque parte do consumo produtivo não
ocorre. A categoria *férias/fim de ano* reúne o período de baixa atividade em torno
da virada do ano, quando muitas indústrias param e a demanda se aproxima da de um
fim de semana prolongado, ainda que seja um dia de semana no calendário. A divisão
é aplicada de forma determinística — o dia da semana e um calendário de feriados
definem o rótulo, sem escolha manual dia a dia — e esse rótulo é usado de duas
maneiras: como característica na formação dos clusters, para que dias operacionalmente
distintos não sejam misturados, e como interpretação do medoide resultante. A coluna
"composição dominante" da tabela informa a estação e o regime majoritários dentro de
cada cluster; por exemplo, "regime 54%" indica que 54% dos dias do bloco D01 caem no
regime de férias/fim de ano, o rótulo herdado pelo medoide. O motivo de fazer isso é
direto: sem esses marcadores, um agrupamento guiado só pela forma horária poderia
fundir um feriado de baixa carga com um fim de semana de baixa carga e, com isso,
distorcer a frequência com que cada regime realmente ocorre no ano — justamente a
frequência que vira peso e anualiza o modelo.

#text(size: 8.8pt)[
#table(
  columns: (0.60fr, 1.05fr, 0.80fr, 0.95fr, 1.25fr, 2.15fr),
  inset: 4pt,
  align: left,
  [Bloco], [Data], [Peso dias], [Estação], [Tipo do dia], [Composição dominante],
  [D01], [30/01/2025], [71], [verão], [férias/fim de ano], [verão 93%; regime 54%],
  [D02], [25/04/2025], [84], [outono], [útil], [outono 99%; útil 69%],
  [D03], [02/08/2025], [104], [inverno], [fim de semana], [inverno 78%; fim de semana 35%],
  [D04], [10/10/2025], [106], [primavera], [útil], [primavera 73%; útil 84%],
)
]

A figura a seguir não traz ainda um resultado econômico; ela é uma verificação de
sanidade da compressão temporal antes de confiarmos nela no resto do relatório.
Queremos enxergar duas coisas: que os quatro medoides cobrem estações distintas e
não deixam uma parte do ano sem representante, e que os pesos de fato somam 365, sem
o que qualquer total anual estaria errado por construção.

#imgfig(
  [Dias representativos e pesos anualizados.],
  "outputs/figures/rel_dias_representativos.png",
  note: [
    Cada barra é um medoide real e sua altura conta quantos dias do ano aquele
    cluster resume; a cor marca a estação e a anotação, o regime de calendário. Os
    pesos 71, 84, 104 e 106 somam exatamente 365, de modo que receitas e custos
    operacionais já saem anualizados, sem multiplicar de novo por 365. O verão ser o
    menor bloco não significa que a estação tenha apenas 71 dias. O algoritmo agrupa
    perfis operacionais, não fatias rígidas do calendário: o cluster de D01 é 93%
    verão, ou cerca de 66 dos seus 71 dias, e os demais dias de dezembro a fevereiro
    ficaram mais próximos de medoides de outros blocos. Isso é plausível porque carga,
    solar, eólica e regime de calendário também entram na distância; dois dias da
    mesma estação podem ser menos parecidos entre si do que dias de estações vizinhas.
    Portanto, a diferença entre 71 e 106 mede a geometria encontrada pelo agrupamento,
    não a duração nem a importância econômica de cada estação.

    A escolha de quatro blocos também merece explicação. No notebook 07, §3, célula
    `[8]`, função `build_representative_days`, ela obtém o
    maior escore de seleção (0,314) e cobre as quatro estações, mas a silhueta de 0,235
    indica separação apenas moderada e os medoides cobrem três de quatro regimes de
    calendário: nenhum feriado aparece como representante. Aumentar para até 12
    clusters tampouco garante capturar esse evento raro, de modo que feriados e
    extremos deveriam entrar em sensibilidades próprias. Há ainda perda de variação
    dentro de cada bloco — ondas de calor, frentes frias e dias nublados atípicos — e
    de sequência entre dias: não há memória de uma seca que se agrava por semanas nem
    de um reservatório enchendo. Por isso cada reservatório é fechado ciclicamente no
    bloco; a escolha evita ganhar energia de graça, mas impede representar transferência
    sazonal de água — uma limitação real, não um detalhe técnico.
  ]
)

O peso mostra quantos dias cada bloco resume, mas não diz se esses dias são realmente
parecidos. Para auditar essa perda de informação, calculamos, para cada membro, a
distância RMS ao centroide no mesmo espaço padronizado que entrou no $k$-médias. Uma
unidade, portanto, significa aproximadamente um desvio-padrão por característica.

#text(size: 8.6pt)[
#table(
  columns: (0.70fr, 0.90fr, 0.90fr, 0.90fr, 0.85fr),
  inset: 4pt,
  align: right,
  [Bloco], [RMS média], [P90 RMS], [Máxima RMS], [Peso],
  [D01], [0,899], [1,252], [1,544], [71],
  [D02], [0,775], [1,015], [1,323], [84],
  [D03], [0,799], [1,016], [1,211], [104],
  [D04], [0,697], [1,052], [1,544], [106],
)
]

#imgfig(
  [Homogeneidade interna dos quatro clusters.],
  "outputs/figures/rel_homogeneidade_clusters.png",
  note: [
    Cada caixa resume quão longe os dias de um cluster ficam de seu medoide: a linha
    central é a mediana, a caixa contém os 50% centrais, os bigodes cobrem a faixa
    usual de 1,5 intervalos interquartis e os pontos isolados são dias especialmente
    distantes. D04 chama atenção por combinar a menor distância média, 0,697, com o
    maior peso, 106 dias. Isso indica que ele representa com boa proximidade a parte
    central de um bloco grande. Seu máximo de 1,544, porém, mostra que há poucos dias
    atípicos mesmo nesse cluster coeso. D01 é o bloco mais heterogêneo: sua média é
    0,899 e seu P90 chega a 1,252, de modo que a cauda superior contém dias menos
    parecidos com o representante. D02 e D03 ficam no meio, com P90 quase idêntico,
    em torno de 1,02.

    O conjunto não sugere clusters completamente homogêneos nem uma compressão
    inútil. Ele sugere uma representação razoável do centro de cada regime, com perda
    relevante nas caudas. Isso é coerente com a silhueta de 0,235: os quatro grupos
    existem, mas se sobrepõem moderadamente. Assim, os medoides são adequados para o
    cenário central anualizado, enquanto feriados, ondas de calor, frentes frias e
    dias renováveis extremos ainda pedem sensibilidades próprias.
    O resultado vem de `representative_days_homogeneity.csv`, produzido no notebook
    07, §3, célula `[8]`, pela função `build_representative_days`; procure o bloco
    que constrói `rep_day_homogeneity`. A figura é reproduzida no notebook 08, §11,
    célula `[12]`. Assim, resultado, implementação e artefato tabular
    apontam para a mesma auditoria.
  ]
)

== Parâmetros, natureza e status de validação

Nem todo número abaixo tem o mesmo estatuto. “Observado” significa lido de uma
fonte externa; “proxy construída” é uma transformação transparente dos dados;
“regulatório” é um parâmetro oficial usado por aproximação; “calibrado” foi escolhido
para reproduzir um alvo; “hipótese estilizada” organiza o laboratório sem pretensão
de estimação; e “regularizador numérico” existe principalmente para estabilidade.
Essa taxonomia impede tratar uma escolha computacional como evidência empírica.

#text(size: 7.35pt)[
#table(
  columns: (1.15fr, 1.45fr, 1.05fr, 2.00fr, 1.20fr),
  inset: 3.0pt,
  align: left,
  [Objeto], [Baseline], [Natureza], [Origem ou justificativa], [Sensibilidade/validação],
  [Perfis $L,G^s,G^w,G^h,G^"th"$ e CMO], [Painel horário de 2025], [Observado],
  [Arquivos do ONS; agregações e lacuna de CMO descritas acima.], [Consistência de chaves e 35.040 linhas verificada.],
  [Capacidades por fonte], [`p99` horário; solar total 39,069 GW], [Proxy construída],
  [Percentil 99 da geração agregada por subsistema; não é placa instalada.], [Substituição por cadastro/disponibilidade pendente.],
  [Tipos solares $theta$], [$xi_"high"=1,08$, custo $0,94q$ e massa 45%; $xi_"low"=0,94$, custo $1,08q$ e massa 55%], [Hipótese estilizada],
  [Qualidade composta e persistente de projeto/desenvolvedor; não estimada.], [Massas, amplitudes e tipos cruzados pendentes.],
  [Custo solar $q_ell$], [N: 514,9; NE: 437,3; SE: 446,0; S: 475,3 mil R\$/MW/ano], [Proxy construída],
  [$("CAPEX"+"conexão") "FRC"+"O&M"$, com WACC 8% e vida de 25 anos.], [Sensibilidade a CAPEX, WACC e conexão pendente.],
  [Déficit $pi_u$], [8.327,76 R\$/MWh], [Parâmetro regulatório usado como proxy],
  [Custo oficial de déficit de 2025, CCEE Comunicado 999/24; não é VOLL estimado.], [Sensibilidade do planejador pendente.],
  [Água $pi_h$], [250 R\$/MWh], [Hipótese estilizada],
  [Cenário central de custo de oportunidade da água.], [Executada em 80, 250 e 500 R\$/MWh.],
  [Curtailment $pi_s,pi_w$], [0,01 R\$/MWh], [Regularizador numérico],
  [Desempate que prefere entregar energia disponível; sem interpretação tarifária.], [Sensibilidade específica pendente.],
  [Curvaturas do clearing], [$epsilon_n=2 times 10^(-4)$, $gamma_r=2 times 10^(-3)$ e $epsilon_x=10^(-2)$], [Regularizadores numéricos],
  [Convexificação, suavização de rampa e unicidade/estabilidade do despacho.], [Resíduos e variação paramétrica executados; ver sensibilidade quadrática.],
  [Reservatório equivalente], [$overline(S)_ell=8K_ell^h$, estoque inicial 50% e fechamento cíclico], [Hipótese estilizada],
  [Evita mineração do estoque dentro dos medoides; não reproduz cascatas.], [Afluência testada em 70%, 100% e 130%.],
  [Rede reduzida], [8,0 GW I--SE; 8,0 GW I--NE; 4,5 GW I--N; 8,0 GW NE--SE; 6,0 GW SE--S], [Hipótese estilizada],
  [Limites fixos do laboratório, não limites horários estimados.], [Calibração externa pendente.],
  [Dias representativos], [$k=4$; pesos 71, 84, 104 e 106], [Seleção de modelo],
  [Silhueta, cobertura sazonal/calendário e parcimônia entre $k=4,dots,12$.], [Critérios comparados; extremos fora da amostra pendentes.],
  [Escala estática $kappa_u$], [244,183 MW por bin], [Normalização numérica],
  [$"39.069,292" " MW"/160$; converte o índice $k$ em MW.], [Reparametrização, não parâmetro econômico.],
  [Curvatura estática $gamma$ e prêmio $p$], [0,312 milhão R\$/bin² e $p=19,1$ mil R\$/MW/ano], [Calibrado],
  [Bisseção para reproduzir 39,069 GW de solar total.], [Total reproduzido; distribuição regional não validada.],
  [Horizonte e desconto], [$H=6$, $beta=0,95$ por estágio], [Hipótese estilizada],
  [Laboratório abstrato; não corresponde automaticamente a seis anos.], [Sensibilidade temporal pendente.],
  [Grade e população dinâmica], [93 MW por bin; $M=105,025$ agentes equivalentes], [Normalização numérica],
  [$M="39.069,292"/(4 times 93)$ fixa a escala inicial.], [Identidade de escala, não contagem de usinas.],
  [Atritos dinâmicos], [$chi=2.000$ R\$/MW²; entrada 20 milhões; mudança 500 milhões; concentração 5 milhões de R\$], [Hipóteses estilizadas],
  [Impedem ajustes e realocações instantâneos e concentração degenerada.], [Sensibilidades pendentes.],
  [Temperatura e relaxação], [$tau=100$ milhões de R\$; $alpha=0,10$ dinâmico e 0,06 estático], [Regularizadores numéricos],
  [Suavizam política e iteração; não são preferências estimadas.], [Ponto fixo baseline certificado; continuação em $tau$ e $alpha$ pendente.],
)
]

#pagebreak(weak: true)

= Mecanismo empírico: carga residual, preço e captura

== Carga líquida, carga residual e rampa

Antes de qualquer modelo de decisão, precisamos de uma linguagem para descrever o
problema que a solar e a eólica criam ao longo do dia. Essa linguagem são três
objetos derivados da carga, e vale defini-los um a um, porque eles significam coisas
diferentes e são frequentemente confundidos.

O ponto de partida é a carga observada $L_(ell,t)$, em MW: é a potência total que o
subsistema $ell$ precisa consumir na hora $t$, antes de descontar qualquer geração.
A geração solar observada é $G_(ell,t)^s$ e a eólica observada é $G_(ell,t)^w$, ambas
também em MW. A partir delas definimos a *carga líquida de solar*,

$ L_(ell,t)^"net" = L_(ell,t) - G_(ell,t)^s, $

que representa quanto ainda falta atender depois de a solar produzir. Subtraindo
também a eólica, chegamos à *carga residual*,

$ L_(ell,t)^"res" = L_(ell,t) - G_(ell,t)^s - G_(ell,t)^w. $

A carga residual é o objeto central: ela aproxima quanta potência o *restante* do
sistema — hidro, térmica, intercâmbio — precisa fornecer numa dada hora, depois que
as fontes variáveis já entregaram tudo o que podiam de graça. É importante frisar o
que ela *não* é. Não é uma demanda econômica no sentido de responder a preço, e não
embute, por si só, restrições de rede, limite de rampa térmica ou estoque de água;
é apenas uma conta de energia, não uma solução de despacho. Ela diz *quanto* o
sistema precisa produzir, não *se consegue* produzir naquele instante.

O terceiro objeto responde a uma pergunta distinta. Enquanto a carga residual é um
*nível* — quanto o sistema precisa entregar naquela hora —, a *rampa residual*

$ Delta L_(ell,t)^"res" = L_(ell,t+1)^"res" - L_(ell,t)^"res" $

é uma *variação entre horas consecutivas*, em MW por hora: quanto o restante do
sistema precisa subir ou descer de uma hora para a seguinte. A diferença é
econômica, não apenas aritmética. Um nível residual alto exige *capacidade* (haver
MW suficientes instalados); uma rampa residual íngreme exige *flexibilidade* (haver
usinas capazes de variar rápido a produção). A solar mexe nas duas ao mesmo tempo,
mas em direções que não coincidem: ela reduz o nível residual ao meio-dia e, ao se
pôr, cria a rampa mais dura do dia. Guardar essa distinção é o que permite, mais à
frente, entender por que um dia com muita energia renovável ainda pode ter aperto no
entardecer.

A primeira pergunta empírica é se o preço reage à carga residual da forma que a
intuição sugere: quando sobra pouca energia barata para atender, o preço deveria
subir. Queremos testar isso diretamente nos dados de 2025, subsistema por
subsistema, para saber se o mecanismo que vamos modelar tem lastro observado ou se é
só teoria. É um exercício exploratório e descritivo — não uma estimação estrutural —,
mas ele decide se vale a pena endogenizar o preço no clearing.

#imgfig(
  [CMO observado versus carga residual por subsistema.],
  "outputs/figures/rel_cmo_vs_carga_residual.png",
  note: [
    Cada ponto é uma hora de 2025; no eixo horizontal está a carga menos solar menos
    eólica e no vertical o CMO observado, com uma reta de regressão descritiva e a
    correlação de Pearson anotada em cada painel. O sinal esperado aparece em três
    das quatro regiões — 0,581 no Norte, 0,412 no Sudeste/CO e 0,281 no Sul —, mas o
    que mais chama atenção é quão *fraca* é essa relação: ao quadrado, essas
    correlações correspondem a apenas 34%, 17% e 8% da variação linear, respectivamente.
    Mesmo no melhor caso, a carga residual deixa a maior parte do preço sem explicar,
    sinal de que água, déficit e rede pesam tanto quanto o nível residual. O Nordeste é o caso
    instrutivo: a correlação é negativa, -0,122, porque a eólica local muitas vezes
    excede a carga regional e o preço passa a ser ditado por excedente e intercâmbio,
    não pela carga própria. Isso é um alerta honesto: a leitura local ingênua
    simplesmente quebra onde a renovável é dominante. Vale sublinhar o limite do
    exercício — correlação não identifica a curva de oferta nem prova causalidade; o
    valor da figura é apenas documentar o canal que o clearing, adiante, tornará
    estrutural em vez de imposto. A correlação negativa do Nordeste tampouco refuta a
    canibalização: mostra que uma regressão contemporânea e local não separa o efeito
    da renovável dos estados hidrológicos e das trocas de energia que ocorrem junto
    com ela.
  ]
)

Também é importante saber se o preço
tem uma forma horária estável, que justificaria trabalhar com uma curva média por
hora, ou se ele é tão disperso dentro de cada hora que uma média esconde mais do que
revela. A resposta orienta uma decisão metodológica concreta — se a dispersão for
grande, não podemos representar o ano por um único perfil de preço, e a compressão
por medoides sazonais passa a ser necessária, não conveniente.

#imgfig(
  [Distribuição horária do CMO observado por subsistema.],
  "outputs/figures/rel_cmo_perfil_horario_boxplot.png",
  note: [
    A linha central de cada caixa é a mediana; a caixa cobre o intervalo interquartil
    e os bigodes se estendem até 1,5 vezes esse intervalo. Os outliers são omitidos
    apenas da figura para preservar a escala, não do cálculo. A dispersão e a
    assimetria dentro de uma mesma hora confirmam que uma curva média única esconderia
    regimes muito diferentes. A figura é gerada no notebook 08, §2, célula `[3]`,
    a partir do painel preparado no notebook 07, §1.
  ]
)

A curva do pato torna nível e rampa residual visíveis de uma vez só. Os quatro
subsistemas aparecem lado a lado porque não basta ter energia: é preciso tê-la na
hora e no nó corretos.

#imgfig(
  [Curvas do pato por subsistema: carga, líquida de solar e residual.],
  "outputs/figures/rel_curva_pato.png",
  note: [
    A curva preta é a carga média; a laranja subtrai a solar; a verde subtrai também
    a eólica. Os vales residuais médios ocorrem às 12h no SE (30,06 GW) e no S
    (9,23 GW), às 7h no NE (-5,19 GW) e às 12h no N (6,17 GW). As maiores rampas
    horárias são 5,17 GW às 17h no SE, 1,32 GW às 16h no S, 1,78 GW às 17h no NE
    e 0,49 GW às 14h no N. O residual negativo no NE significa excedente médio local
    de solar e eólica e, portanto, necessidade de exportação ou corte. São perfis
    médios; dias individuais podem ser mais severos. A figura e esses cálculos estão
    no notebook 08, §3, célula `[4]`.
  ]
)

== Receita capturada e canibalização

A seção anterior estabeleceu que o preço cai quando a solar rebaixa a carga
residual. Falta juntar as duas pontas: *de quem* é a receita afetada por isso? A
resposta é o próprio gerador solar. Um projeto não vende energia de forma uniforme
ao longo do dia; ele vende quando produz, isto é, exatamente nas horas de sol. Se
são justamente essas as horas em que a solar, no agregado, derruba o preço, então o
gerador colhe um preço médio pior do que a média do sistema — e quanto mais solar
entra, pior fica. Esse é o fenômeno de canibalização de receita, e o objetivo desta
seção é medi-lo nos dados e mostrar que a matemática do modelo o reproduz.

Para quantificá-lo precisamos de uma medida de "preço que a fonte de fato recebe".
Seja $G_t^r$ a geração da fonte $r$ (solar ou eólica) na hora $t$, em MWh naquela
hora, e $P_t$ a série de preço, em R\$/MWh. O *preço capturado* é o preço médio
ponderado pela própria geração,

$ P^"cap,r" = frac(sum_t P_t G_t^r, sum_t G_t^r), $

também em R\$/MWh. Ele responde à pergunta: por cada MWh gerado ao longo do ano,
quanto valeu, em média, a hora em que ele foi gerado? Como só há peso nas horas em
que a fonte produz, uma fonte que gera em horas caras captura mais do que a média e
uma que gera em horas baratas captura menos. Para comparar fontes e regiões numa
escala única, dividimos o preço capturado pela média simples do preço no período,
$T^(-1) sum_t P_t$, obtendo o *fator de captura* adimensional

$ F^"cap,r" = frac(P^"cap,r", T^(-1) sum_t P_t). $

A leitura é direta: $F^"cap,r"=1$ significa que a fonte recebe exatamente a média do
sistema; abaixo de um significa que ela concentra geração em horas de preço abaixo
da média — o desconto de canibalização; acima de um, o contrário. No diagnóstico
empírico usamos $P_t = CMO_t^"obs"$, o custo marginal de operação observado; convém
repetir que isso é um termômetro de custo sistêmico, não o PLD que remunera
contratos.

O mesmo fenômeno pode ser derivado, e não apenas medido. Escreva
$lambda_t = Lambda(L_t - K^"agg" a_t^s - G_t^w)$. $Lambda$ não é uma função
analítica imposta nem uma regressão do CMO: é a forma reduzida do mapa implícito
produzido pelo clearing convexo,
$Lambda: L mapsto partial V^*(L)/partial L$. Como o valor ótimo $V^*$ é convexo na
carga, esse mapa é monotônico em cada direção; no programa quadrático, ele é
suave por trechos e muda de inclinação quando blocos, linhas ou bandas se tornam
ativos. Assim, $Lambda'_t$ é uma derivada marginal local do *modelo*, não uma
constante estrutural e não o coeficiente de uma correlação observada. O termo
$K^"agg" a_t^s$ é a geração solar agregada: $K^"agg"$ é a capacidade solar total
instalada no sistema, em MW, e $a_t^s$ é o *perfil solar normalizado*, a fração
adimensional dessa capacidade efetivamente disponível na hora $t$ (perto de zero à
noite, próximo do máximo ao meio-dia). O produto dos dois, portanto, é a potência
solar em MW.

Um projeto pequeno tem seu próprio perfil $a_(i,t)^s$ e recebe o preço capturado

$ R_i^"cap" = frac(sum_t lambda_t a_(i,t)^s, sum_t a_(i,t)^s). $

Mantendo os perfis fixos e derivando em relação à capacidade solar agregada do
sistema,

$ frac(partial R_i^"cap", partial K^"agg")
  = - frac(sum_t a_(i,t)^s a_t^s Lambda'_t, sum_t a_(i,t)^s) <= 0. $

O sinal é inequivocamente não positivo: adicionar capacidade solar ao sistema nunca
aumenta, e em geral reduz, o preço capturado de um projeto solar já existente. A
desigualdade é estrita — a canibalização morde de verdade — quando o projeto produz
nas mesmas horas em que a solar agregada produz ($a_(i,t)^s a_t^s > 0$) e nessas
horas o preço responde à carga residual ($Lambda'_t > 0$). Nada disso depende de o
projeto ser grande: cada gerador é tomador de preço, mas todos juntos movem o preço
que cada um recebe — a tensão exata que o jogo de campo médio vai formalizar.

No modelo, $Lambda$ não é uma função escolhida por regressão. É o mapa implícito que
leva a carga residual ao dual do balanço de energia depois que o QP é resolvido.
Assim, $Lambda'$ representa quanto o preço-sombra marginal se inclina quando a carga
residual aumenta, mantendo os demais dados fixos. Em um trecho suave, ele pode ser
lido como derivada local; em uma quina, quando muda o conjunto de restrições ativas,
a inclinação pode saltar.

No código, usamos uma diferença finita direcional. Depois de resolver o caso de
referência, acrescentamos 100 MW à carga de cada região, uma por vez, reexecutamos
exatamente o mesmo QP e calculamos

$ hat(Lambda)'_ell(Delta L) =
  frac(overline(lambda)_ell(L+Delta L e_ell)-overline(lambda)_ell(L),Delta L),
  quad Delta L=100 " MW". $

O objeto calculado é, portanto, uma *secante direcional de 100 MW*. Ele aproxima uma
derivada apenas quando esse intervalo permanece no mesmo regime marginal.

#text(size: 8.5pt)[
#table(
  columns: (0.75fr, 1.55fr),
  inset: 4pt,
  align: right,
  [Subs.], [$hat(Lambda)'(100)$ R\$/MWh por MW],
  [SE], [0,000286],
  [S], [0,187593],
  [NE], [0,000293],
  [N], [0,000293],
)
]

#imgfig(
  [Sensibilidade secante do preço-sombra a mais 100 MW de carga.],
  "outputs/figures/rel_lambda_prime_sensitivity.png",
  note: [
    O Sul chama atenção porque 0,187593 é cerca de 650 vezes o valor das outras
    regiões. Esse número não significa que todo MW adicional no Sul sempre aumente o
    preço nessa proporção. A perturbação de 100 MW atravessa uma quina do despacho.
    No medoide D02, que representa 84 dias do ano, a interface SE--S já está no limite
    em 21 das 24 horas no caso de referência; com a perturbação, ela fica saturada nas
    24 horas. Ao mesmo tempo, o primeiro bloco térmico sulista está no limite e passa
    a ser necessário acionar o bloco seguinte, de custo marginal maior. O preço-sombra
    médio do Sul em D02 sobe 81,40 R\$/MWh, enquanto as variações nos demais medoides
    ficam entre 0,019 e 0,048 R\$/MWh. Como D02 recebe peso de 84 dias, ele domina a
    média anual.

    A auditoria por tamanho de passo confirma a quina: para perturbações de 1 a 10 MW,
    a secante sulista fica perto de 0,0011 R\$/MWh por MW; entre 60 e 65 MW, a última
    folga da interface desaparece e o valor salta. Logo, 0,187593 descreve a resposta
    média do intervalo de 100 MW que cruza dois regimes, não a inclinação infinitesimal
    no ponto de referência. SE, NE e N não atravessam mudança comparável de rede e
    mérito térmico nesse mesmo teste, por isso seus valores permanecem próximos.

    A equação principal é implementada no notebook 07, §4.6, célula `[14]`: procure
    `lambda_prime_rows`, `perturbed_load[:, li, :] += delta_load_MW` e
    `DispatchEngineCVX(...).solve`. A auditoria do Sul está no mesmo bloco, em
    `south_lambda_step_rows`, e gera `lambda_prime_south_step_sensitivity.csv`.
    O resultado regional é salvo em `lambda_prime_sensitivity.csv` e a figura é
    reproduzida no notebook 08, §15, célula `[16]`.
  ]
)

#text(size: 8.8pt)[
#table(
  columns: (0.72fr, 1.15fr, 1.15fr),
  inset: 4pt,
  align: right,
  [Subs.], [Fator solar], [Fator eólico],
  [SE], [0,725], [1,044],
  [S], [0,729], [1,067],
  [NE], [0,640], [1,165],
  [N], [0,684], [1,309],
)
]

A derivação diz que a solar *deveria* capturar abaixo da média; a figura seguinte
verifica se isso já aparece nos dados de 2025, antes mesmo de qualquer expansão
adicional, e se a eólica se comporta de forma distinta. Se o desconto solar já
existe hoje, ele é um ponto de partida, não uma previsão — e isso muda como se deve
ler qualquer projeção de receita por MW.

#imgfig(
  [Fatores de captura solar e eólico usando CMO observado.],
  "outputs/figures/rel_fatores_captura.png",
  note: [
    A linha de referência em um separa quem captura acima da média local do CMO de quem
    captura abaixo. O resultado é nítido e na direção prevista: a solar fica abaixo
    de um nos quatro subsistemas, entre 0,640 e 0,729, enquanto a eólica fica acima,
    entre 1,044 e 1,309. Em números, a solar do Nordeste recebe cerca de 36% menos do
    que a média local do CMO — um desconto grande, que já existe *sem* a expansão que o
    modelo vai estudar, sinal de que o ponto de partida já é desfavorável e tende a
    piorar. Vale notar o contraste locacional: o desconto solar é mais fundo no
    Nordeste (0,640) do que no Sudeste/CO (0,725), o oposto do que a distribuição de
    capacidade sugeriria, e a eólica nordestina é a que mais captura acima da média —
    a complementaridade temporal entre vento e carga trabalha a favor dela na
    amostra. Isso não prova que a eólica esteja imune à canibalização: prova apenas
    que, com a capacidade e os perfis de 2025, ela produziu em horas cujo CMO ficou
    acima da média local; uma expansão eólica pode deslocar justamente essas horas e
    reduzir o fator. A cautela central é que isto é um fator calculado com CMO, não um retorno
    financeiro realizado: contratos de longo prazo, PLD, encargos e a diferença entre
    geração centralizada e distribuída podem amortecer boa parte desse desconto na
    receita efetiva de um projeto. O número mede a pressão de mercado, não o resultado
    de caixa.
  ]
)

#resultbox([Fato estilizado que alimenta o modelo])[ 
  O dado e a derivação apontam na mesma direção: ao crescer, a solar reduz a carga
  residual nas próprias horas de produção e tende a reduzir seu preço capturado.
  O restante do relatório endogeniza essa relação por meio do clearing, em vez de
  impor uma curva de preço arbitrária.
]

#pagebreak(weak: true)

= Clearing locacional compartilhado

Toda a parte empírica descreveu o problema; agora precisamos de um mecanismo que o
*resolva*. Esse mecanismo é o clearing: dado o parque instalado, ele decide, hora a
hora e região a região, quanto cada usina gera para atender a carga ao menor custo,
respeitando rede, água e limites térmicos, e devolve o preço-sombra que sustenta
essa decisão. É o coração físico do relatório, e a escolha metodológica mais
importante é que *todas* as camadas de decisão adiante — o laboratório de dois
produtores, o benchmark privado, o contrafactual, o planejador social e o MFG —
chamam exatamente este mesmo operador.

Formalmente, o clearing é uma transformação que parte de um vetor de capacidade
solar $K$ (por localização e tipo, em MW) e produz todas as variáveis de operação:

$ K arrow.r overline(g)^s(K) arrow.r
  (g^s,c^s,g^w,c^w,h,S,v,n,F,u,lambda). $

Lida da esquerda para a direita: a capacidade $K$ determina primeiro quanta solar
está *disponível*, $overline(g)^s(K)$; e dessa disponibilidade, junto com as demais
tecnologias, alimenta o problema de despacho, cuja solução é a tupla à direita —
solar entregue e cortada, eólica entregue e cortada, hidro, estoque e vertimento,
geração térmica, fluxos na rede, déficit e, por fim, o preço $lambda$. Cada um
desses símbolos é definido nas subseções a seguir, na ordem em que aparece no
problema.

O clearing é, matematicamente, um *programa quadrático convexo* (QP). “Primal” é o
problema escrito nas quantidades físicas que o operador escolhe — geração, corte,
estoque, vertimento, fluxo e déficit. Sua função objetivo soma termos lineares e
quadráticos convexos, enquanto balanços e dinâmicas são igualdades afins e os
limites são desigualdades lineares. “Dual” é o problema associado que atribui um
multiplicador a cada restrição; esses multiplicadores medem quanto o valor ótimo
muda quando se relaxa marginalmente o respectivo limite. O $lambda$ do balanço
nodal é um desses multiplicadores, depois da normalização pelo peso do medoide.

As condições de Karush--Kuhn--Tucker (KKT) são o sistema que combina quatro
exigências: viabilidade primal (as quantidades obedecem às restrições), viabilidade
dual (multiplicadores de desigualdades têm o sinal correto), estacionariedade (não
há direção marginal de redução do custo depois de incorporar os multiplicadores) e
complementaridade (uma desigualdade folgada tem multiplicador zero; um multiplicador
positivo só pode acompanhar uma restrição ativa). Como o clearing é convexo, com
função diferenciável e restrições afins/caixas factíveis, as KKT são necessárias e
suficientes para caracterizar seu ótimo sob a regularidade usual. É por isso que
podemos usar as condições abaixo para interpretar economicamente $lambda$; essa
leitura não seria automática em um problema não convexo.

== Rede e variáveis de decisão

Há quatro nós físicos e um hub contábil $I$. As arestas são

$ I arrow.l.r "SE", quad I arrow.l.r "NE", quad I arrow.l.r "N",
  quad "NE" arrow.l.r "SE", quad "SE" arrow.l.r "S". $

Para transformar cada fluxo orientado em importação líquida, o código fixa a ordem
e a orientação positiva

$ e_1=(I,"SE"), quad e_2=(I,"NE"), quad e_3=(I,"N"), quad
  e_4=("NE","SE"), quad e_5=("SE","S"). $

Assim, $F_(d,e,t)>0$ significa energia fluindo da primeira para a segunda ponta do
par; um valor negativo representa o sentido inverso. Definindo a matriz de
incidência nó--aresta por

$
  A_(ell e)=cases(
    +1 & "se a aresta " e " entra em " ell,
    -1 & "se a aresta " e " sai de " ell,
    0 & "caso contrário",
  ),
$

o fluxo líquido de importação que aparece no balanço é

$ F_(d,ell,t)^"net"=sum_(e in cal(E)) A_(ell e)F_(d,e,t). $

Com as orientações acima, isso equivale a

$ F_"SE"^"net"=F_(e_1)+F_(e_4)-F_(e_5), quad
  F_"NE"^"net"=F_(e_2)-F_(e_4), quad
  F_"N"^"net"=F_(e_3), quad
  F_"S"^"net"=F_(e_5), $

com os índices $(d,t)$ omitidos apenas para aliviar a notação. No hub,
$F_I^"net"=-F_(e_1)-F_(e_2)-F_(e_3)=0$. Essa última igualdade é a conservação no
nó contábil; ela impede que o hub crie ou consuma energia.

O hub não possui carga nem geração e satisfaz conservação algébrica dos fluxos. Os
limites reduzidos são 8,0 GW em I-SE, 8,0 GW em I-NE, 4,5 GW em I-N,
8,0 GW em NE-SE e 6,0 GW em SE-S. Não existem ligações diretas N-NE ou N-SE no
grafo reduzido.

No código, a definição é única e comum às duas contas: `EDGES` contém
`("NE","SE")`, `FALLBACK_FMAX` fixa seus 8.000 MW e `N_EDGES` é derivado do
comprimento dessa lista. Tanto `DispatchEngineCVX._build_problem` quanto
`solve_social_planner` percorrem os mesmos pares; portanto, balanço, limites, clearing
privado e planejador usam a ligação direta sem duplicar uma equação paralela.

#nbref([
  Notebook 07, §4.1, célula `[10]`: procure `EDGES`, `FALLBACK_FMAX`,
  `N_EDGES` e `FMAX`. O balanço operacional que percorre essas arestas está em
  §4.2, célula `[11]`, classe `DispatchEngineCVX`, método `_build_problem`. A
  mesma lista entra no planejador em §8, célula `[22]`, função
  `solve_social_planner`, nos laços
  `for edge, (origin, destination) in enumerate(EDGES)`.
])

Antes de escrever as equações, vale visualizar a rede sobre a qual elas operam,
porque é ela que transforma o problema de "quanto gerar" em um problema *locacional*:
energia barata no Nordeste só ajuda o Sudeste/CO se houver linha e folga para
transportá-la. A figura mostra a topologia reduzida que o clearing enxerga, e o
objetivo de olhá-la é fixar quais trocas são possíveis e onde o congestionamento
pode isolar uma região.

#imgfig(
  [Topologia reduzida de intercâmbio utilizada pelo clearing.],
  "outputs/figures/rel_topologia_intercambio.png",
  note: [
    Os quatro círculos coloridos são subsistemas com carga e geração próprias; o
    círculo cinza $I$ é um nó puramente contábil, sem carga nem geração, que apenas
    conserva os fluxos. Cada seta dupla é *uma* variável de fluxo com sinal e um
    limite bilateral — não dois fluxos independentes. Nordeste e Sudeste/CO possuem
    uma ligação direta de 8,0 GW; o Norte continua ligado apenas pelo hub e o Sul,
    pelo Sudeste/CO. Essa é uma
    representação reduzida: cinco
    interfaces com limites fixos capturam o congestionamento inter-regional de
    primeira ordem, mas ignoram perdas, restrições internas de cada subsistema,
    contingências e a variação horária dos limites que o ONS de fato pratica.
    Conclusões locacionais herdam essa simplificação.
  ]
)

Para cada $(d,ell,t)$, as variáveis do primal são:

- $g^s$ e $c^s$: solar entregue e solar cortada;
- $g^w$ e $c^w$: eólica entregue e eólica cortada;
- $h$, $S$ e $v$: geração hídrica, estoque e vertimento;
- $n_b$: geração do bloco térmico $b$;
- $u$: déficit, isto é, energia não atendida penalizada;
- $F_e$: fluxo orientado na aresta $e$.

== Pilha térmica e geração nuclear exógena

Cada subsistema possui exatamente três blocos térmicos, baixo, médio e alto. A
construção parte da mediana regional do CMO horário já preenchido,

$ c_ell^"base"="clip"("median"_t(CMO_(ell,t)^"obs"),60,450), $

e aplica multiplicadores de custo $(0,80;1,25;1,85)$. A potência térmica `p99`
$K_ell^"th,p99"$ é repartida nas frações $(0,50;0,30;0,20)$. Portanto,

$ (c_(ell,1)^"th",c_(ell,2)^"th",c_(ell,3)^"th")
  =c_ell^"base"(0,80;1,25;1,85), $

$ (overline(N)_(ell,1),overline(N)_(ell,2),overline(N)_(ell,3))
  =K_ell^"th,p99"(0,50;0,30;0,20). $

As quatro medianas observadas ficam dentro do intervalo de truncamento; portanto,
o `clip` não altera nenhuma delas nesta execução. Os valores efetivamente enviados
ao clearing são:

#text(size: 8.1pt)[
#table(
  columns: (0.64fr, 0.58fr, 1.00fr, 1.15fr, 1.18fr),
  inset: 3.4pt,
  align: right,
  [Subs.], [Bloco], [$c_ell^"base"$ R\$/MWh], [$overline(N)_(ell,b)$ MW], [$c_(ell,b)^"th"$ R\$/MWh],
  [SE], [1], [263,105], [3.975,254], [210,484],
  [SE], [2], [263,105], [2.385,152], [328,881],
  [SE], [3], [263,105], [1.590,101], [486,744],
  [S], [1], [272,338], [1.025,362], [217,870],
  [S], [2], [272,338], [615,217], [340,422],
  [S], [3], [272,338], [410,145], [503,824],
  [NE], [1], [130,923], [1.484,962], [104,738],
  [NE], [2], [130,923], [890,977], [163,653],
  [NE], [3], [130,923], [593,985], [242,207],
  [N], [1], [175,185], [1.720,439], [140,148],
  [N], [2], [175,185], [1.032,264], [218,981],
  [N], [3], [175,185], [688,176], [324,092],
)
]

Esses custos não são ofertas observadas das usinas. São proxies construídas a partir
de uma única estatística do CMO por subsistema — a mediana de 2025 — e dos
multiplicadores declarados acima.

Já $g_(d,ell,t)^"nuc"$ é a geração nuclear exógena observada no medoide. Ela entra
no balanço como um perfil fixo `nuc_avail`, obtido da coluna nuclear do painel, e não
é variável do primal: o clearing não decide seu despacho nem seu curtailment. No
painel utilizado, a geração nuclear aparece apenas no Sudeste/CO.

== Oferta renovável e curtailment explícito

A primeira etapa da transformação é converter capacidade em energia disponível. Para
a solar, a potência que *poderia* ser gerada no medoide $d$, na localização $ell$ e
na hora $t$ é

$ overline(g)_(d,ell,t)^s(K)
  = a_(d,ell,t)^s sum_theta xi_theta K_(ell,theta). $

Cada símbolo tem um papel concreto. A barra sobre $g$ indica que $overline(g)^s$ é um
*teto* de geração, não a geração efetiva — o máximo que o sol permite naquela hora,
em MW. O perfil $a_(d,ell,t)^s$ é a fração adimensional da capacidade disponível
(zero à noite, próximo de um ao meio-dia), o mesmo objeto da seção empírica. A soma
$sum_theta xi_theta K_(ell,theta)$ agrega os dois tipos solares da região: $K_(ell,
theta)$ é a capacidade instalada do tipo $theta$ em MW e $xi_theta$ é o fator de
produtividade daquele tipo — $xi_"high"=1,08$ e $xi_"low"=0,94$ —, de modo que um MW
do tipo `high` rende mais energia disponível do que um MW do tipo `low`. É exatamente
por aqui que os tipos solares definidos na Seção 2 entram no despacho.

Tendo o teto, o despacho decide quanto dessa energia é de fato entregue e quanto é
descartado, separando as duas parcelas de forma explícita:

$ g_(d,ell,t)^s + c_(d,ell,t)^s = overline(g)_(d,ell,t)^s(K), $

$ g_(d,ell,t)^w + c_(d,ell,t)^w = overline(g)_(d,ell,t)^w, $

onde $g^s$ (entregue) e $c^s$ (cortada, o *curtailment*) são ambas variáveis não
negativas em MW, e o análogo vale para a eólica. Tratar o corte como variável própria
— e não como uma subtração calculada no fim — é uma escolha que importa: garante que
o curtailment nunca fique negativo por acidente numérico e faz com que ele entre nas
condições de complementaridade, ou seja, o modelo só corta renovável quando há razão
econômica (excesso local sem escoamento), e o preço-sombra registra isso.

== Da hidráulica física ao reservatório em energia equivalente

A hidro é a peça mais delicada do sistema brasileiro, e vale ser explícito sobre o
que estamos simplificando. Na realidade física, a potência de uma usina $j$ é

$ P_(j,t)^h = eta_(j,t) rho^w g H_(j,t)^"net" Q_(j,t)^"tur", $

isto é, o produto do rendimento do conjunto turbina-gerador $eta$, da densidade da
água $rho^w$, da gravidade $g$, da *altura de queda líquida* $H^"net"$ e da vazão
turbinada $Q^"tur"$. A altura líquida, por sua vez,

$ H_(j,t)^"net" = Z_j^"mon"(V_(j,t))
  - Z_j^"jus"(Q_(j,t)^"tur"+Q_(j,t)^"spill")
  - H_(j,t)^"loss", $

depende do nível de montante $Z^"mon"$ (que sobe com o volume $V$), do nível de
jusante $Z^"jus"$ (que sobe com a vazão defluente) e das perdas hidráulicas
$H^"loss"$. O ponto de trazer essas fórmulas não é usá-las — é deixar claro *o
quanto* estamos deixando de fora: reproduzi-las por usina exigiria curvas
cota-volume-vazão, rendimentos, defluências mínimas, cascatas e restrições
ambientais, dados que não temos nesta escala.

Por isso o modelo colapsa toda a hidráulica de um subsistema em um único
*reservatório equivalente*, contabilizado diretamente em energia. Nessa
representação, $S_(d,ell,t)$ é o estoque de água expresso em MWh-equivalentes (a
energia que ainda pode ser gerada com a água armazenada), $A_(d,ell,t)$ é a afluência
que chega, também já convertida em MWh-equivalentes, $h_(d,ell,t)$ é a geração
hídrica efetiva (o turbinamento, em MWh na hora) e $v_(d,ell,t)$ é o vertimento —
água que passa sem gerar. Como $h = P^h Delta t$, a dinâmica do estoque vira uma
simples conta de "entra menos sai":

No baseline, $A$ não é a ENA observada nem a conversão direta das vazões do arquivo
hidrológico. Para cada medoide e subsistema, primeiro somamos as 24 horas de geração
hídrica observada, limitamos o orçamento a $24K_ell^"h,p99"$ e dividimos esse total
uniformemente entre as 24 horas. Assim, a energia hídrica diária do dia observado
ancora o volume afluente, mas o perfil intradiário é decidido pelo reservatório. A
vantagem é separar “quanta água-equivalente há no dia” de “quando turbinar”; o custo
é apagar a forma horária e a origem física da afluência. Os arquivos diários de
nível, vazão, ENA e EAR servem aqui como referência e rota de extensão, não como
entrada direta desta execução.

$ S_(d,ell,t+1)=S_(d,ell,t)+A_(d,ell,t)-h_(d,ell,t)-v_(d,ell,t). $

O estoque e o turbinamento têm limites, e o dia é fechado ciclicamente:

$ 0 <= S_(d,ell,t) <= overline(S)_ell, quad
  0 <= h_(d,ell,t) <= overline(H)_ell, quad
  v_(d,ell,t) >= 0, $

$ S_(d,ell,24)=S_(d,ell,0)=0,5 overline(S)_ell. $

A capacidade do reservatório equivalente $overline(S)_ell$ é fixada em oito horas da
potência hídrica proxy — uma escolha calibrada, não física. A condição cíclica
$S_24 = S_0$ é o que impede a trapaça de um medoide gerar energia "de graça"
esvaziando um estoque inicial que ninguém precisou encher: ele permite arbitragem
*dentro* do dia (guardar água do meio-dia para o pico da noite), mas proíbe uma
transferência sazonal fictícia entre medoides que, como vimos na Seção 2, seriam
dias desconexos no tempo.

Uma restrição dinâmica escrita no papel só vale se o solver realmente a respeitar.
A figura seguinte não busca um resultado econômico; ela existe para *auditar* que a
equação de movimento e o fechamento cíclico são de fato obedecidos, e para tornar
tangível o que significa "arbitragem intradiária" no reservatório equivalente.

#imgfig(
  [Estado do reservatório, afluência, turbinamento e vertimento.],
  "outputs/figures/rel_reservoir_state.png",
  note: [
    O painel superior segue o estoque $S_t$ do instante 0 ao 24 em um medoide do
    Sudeste/CO; o inferior põe na mesma escala os três fluxos que o movem — afluência
    $A_t$, turbinamento $h_t$ e vertimento $v_t$. Duas coisas merecem leitura. A
    primeira é econômica: o estoque desce quando o turbinamento supera a afluência e
    sobe no caso contrário, e é assim que a hidro desloca energia das horas de sol
    para o pico do entardecer — a "bateria" barata do sistema em ação. A segunda é de
    controle de qualidade e é o motivo real da figura: a curva retorna *exatamente*
    ao estoque inicial e o erro máximo da identidade
    $S_(t+1)-S_t-A_t+h_t+v_t$ é de $5,47 times 10^(-11)$ MWh, ou seja, ruído de
    ponto flutuante. Isso é uma auditoria, não uma ilustração: confirma que o
    fechamento cíclico está ativo e que nenhuma energia aparece do nada. O que a
    figura *não* autoriza é ler $S_t$ como o volume de um reservatório real — é
    energia equivalente agregada de todo o subsistema.
  ]
)

== Problema de despacho

O que o clearing minimiza é o *custo total de operar o sistema* ao longo do ano, e
cada parcela da função objetivo é um custo econômico ou uma penalidade que
representa uma preferência operacional. Escrevendo $p_d$ para o peso do medoide (os
dias que ele representa, o que anualiza o custo),

$
  min sum_d p_d [
    sum_ell sum_(t=0)^23 [
      sum_b (c_(ell,b)^"th" n_(d,ell,b,t)
        + frac(epsilon_n,2) n_(d,ell,b,t)^2)
      + pi_h h_(d,ell,t)
      + pi_s c_(d,ell,t)^s
      + pi_w c_(d,ell,t)^w
      + pi_u u_(d,ell,t)
    ]
    + sum_ell sum_(t=1)^23 gamma_r
      (n_(d,ell,t)-n_(d,ell,t-1))^2
  ] + epsilon_x ||x||_2^2.
$

A geração térmica total que aparece na rampa é, por definição,

$ n_(d,ell,t):=sum_(b=1)^3 n_(d,ell,b,t). $

O somatório de rampa começa em $t=1$ e termina em $t=23$. Portanto, a hora zero não
é comparada com uma geração inicial exógena nem com a hora 23 do mesmo medoide; o
perfil térmico não tem fechamento cíclico. A primeira penalidade é entre as horas
zero e um, e existem 23 diferenças consecutivas por dia. Essa é uma assimetria em
relação ao reservatório, cujo estoque fecha ciclicamente.

Lendo termo a termo: $c_(ell,b)^"th" n_(d,ell,b,t)$ é o custo do combustível
térmico, o preço $c_b^"th"$ (R\$/MWh) do bloco vezes o quanto ele gera $n_b$ (MWh) —
é o termo que ordena o mérito; a pequena parcela $frac(epsilon_n,2) n^2$ apenas
convexifica cada bloco. O termo $pi_h h$ cobra um custo de oportunidade $pi_h$ por
turbinar água agora em vez de guardá-la. As parcelas $pi_s c^s$ e $pi_w c^w$
penalizam levemente o corte solar e eólico, e $pi_u u$ penaliza *pesadamente* o
déficit $u$ — a energia não atendida — com um valor próximo ao custo social de não
atender a carga. O termo $gamma_r (n_t - n_(t-1))^2$, aplicado apenas para
$t=1,dots,23$, desincentiva variações bruscas da térmica entre horas, representando
o custo de rampa. Por fim, $epsilon_x ||x||_2^2$
é um regularizador numérico sobre o vetor $x$ de todas as variáveis de despacho, que
existe para tornar a solução única e estável.

#nbref([
  Notebook 07, §4.2, célula `[11]`, classe `DispatchEngineCVX`. As variáveis são
  criadas em `_build_problem`; procure `objective`, o termo
  `p["thermal_ramp_penalty"] * cp.sum_squares(...)` e
  `cp.Problem(cp.Minimize(objective), constraints)`. As equações e a interpretação KKT aparecem em
  §4.3; a auditoria numérica do reservatório está na célula `[12]`, que grava
  `reservoir_audit.csv`.
])

O parâmetro $pi_u$ representa o *Value of Lost Load* (VOLL), ou valor econômico
atribuído à carga não atendida. No modelo, $u_(d,ell,t)$ é o déficit de potência na
hora $t$; como os intervalos têm duração de uma hora, seu valor em MW corresponde
numericamente à energia não suprida em MWh naquele intervalo. O produto $pi_u u$
é, portanto, o custo que a função objetivo atribui a deixar de atender essa energia.

Este relatório não estima diretamente o VOLL dos consumidores. Utilizamos
$pi_u=8.327,76$ R\$/MWh, o custo oficial de déficit divulgado no
#link("https://www.ccee.org.br/en/web/guest/-/co-divulgacao-do-custo-de-deficit-e-memoria-de-calculo-referente-ao-ano-2025")[Comunicado CCEE 999/24] para 2025, como uma
proxy operacional para o valor da carga perdida. Essa distinção é importante: o
VOLL procura representar o dano econômico da interrupção, enquanto o custo oficial
de déficit é um parâmetro regulatório adotado aqui para penalizar a energia não
suprida. Ele não é o CMO observado nem o PLD.

Quando $u_(d,ell,t)>0$, a condição de complementaridade faz o preço-sombra
$lambda_(d,ell,t)$ aproximar-se de $pi_u$, salvo as pequenas correções introduzidas
pelos regularizadores. Por isso, o VOLL funciona também como um teto econômico
aproximado para o preço marginal do clearing. Quanto maior $pi_u$, mais o modelo
prefere acionar geração cara, importar energia ou investir em capacidade antes de
aceitar déficit.

Os valores do baseline explicitam as prioridades: $pi_u=8.327,76$ R\$/MWh (déficit é
caríssimo), $pi_h=250$ R\$/MWh (água tem preço relevante), $pi_s=pi_w=0,01$ R\$/MWh
(corte quase de graça), $epsilon_n=2 times 10^(-4)$, $gamma_r=2 times 10^(-3)$ e
$epsilon_x=10^(-2)$. É importante ler os três últimos com ceticismo: eles são
*epsilons de desempate*, escolhidos por estabilidade numérica e não por medição, e
não representam compensação regulatória de curtailment nem um custo de rampa
observado. Como o regularizador entra no preço marginal, qualquer conclusão
normativa precisa de sensibilidade a esses valores — ponto retomado nas limitações.

Tudo isso está sujeito ao *balanço nodal*, a restrição de que oferta e demanda se
igualam em cada nó e hora. Escrito na convenção "carga menos geração igual a zero",

$
  L_(d,ell,t)
  - g_(d,ell,t)^s - g_(d,ell,t)^w - h_(d,ell,t)
  - sum_b n_(d,ell,b,t) - g_(d,ell,t)^"nuc"
  - F_(d,ell,t)^"net" - u_(d,ell,t) = 0,
$

ele soma toda a geração local (solar, eólica, hidro, térmica somada nos blocos e
nuclear), o fluxo líquido de importação $F^"net"$ e o déficit $u$, e exige que
cubram a carga $L$. A escolha do sinal aqui não é cosmética: como veremos, é ela que
faz o preço-sombra sair com o sinal certo, sem correção posterior. Os blocos térmicos
respeitam $0<=n_(ell,b,t)<=overline(N)_(ell,b)$, os fluxos respeitam
$-overline(F)_e<=F_(e,t)<=overline(F)_e$ e o hub contábil fecha com soma algébrica
nula.

== Dual, teorema do envelope e KKT

O clearing produz um preço, mas um número só é útil se soubermos o que ele significa.
O objetivo desta subseção é abrir o preço-sombra e mostrar que ele não é uma caixa
preta: é a soma legível de mérito térmico, valor da água, penalidade de déficit e
congestionamento. A ferramenta para isso são as condições de otimalidade (KKT) e o
teorema do envelope, que ligam o multiplicador do balanço ao custo de uma unidade a
mais de carga. Fazemos isso porque é esse preço decomposto que, mais adiante,
remunera o investidor e conecta a operação à decisão de entrar ou não no mercado.

O multiplicador do balanço nodal é o preço-sombra. Chamando $tilde(lambda)_(d,ell,t)$
o dual bruto dessa restrição e lembrando que o objetivo do medoide foi ponderado por
$p_d$, o preço reportado é normalizado desfazendo esse peso,

$ lambda_(d,ell,t)=frac(tilde(lambda)_(d,ell,t),p_d). $

Com o balanço escrito como carga menos geração igual a zero, o teorema do envelope
fornece diretamente

$ lambda_(d,ell,t)=frac(partial V^*, partial L_(d,ell,t)), $

onde $V^*$ é o valor ótimo do despacho, isto é, o próprio custo total mínimo que a
função objetivo atinge. Em palavras, $lambda$ é exatamente quanto o custo do sistema
sobe se a carga naquele nó e hora aumentar em uma unidade — a definição econômica de
custo marginal. Não há troca posterior de sinal. Em um bloco térmico interior, numa
hora também interior $1<=t<=22$, e omitindo apenas os termos de banda, a condição de
primeira ordem é

$ lambda_t = c_b^"th" + epsilon_n n_(b,t)
  + 2 gamma_r(2n_t-n_(t-1)-n_(t+1))
  + "termo do regularizador". $

Nas bordas, a contribuição da rampa à estacionariedade é somente
$2gamma_r(n_0-n_1)$ em $t=0$ e $2gamma_r(n_23-n_22)$ em $t=23$. Não aparece
$n_(-1)$ nem $n_24$: essa KKT confirma a mesma condição de borda não cíclica usada
no código.

Para o déficit, complementaridade implica que $u_t>0$ aproxima o preço marginal
de $pi_u$, com a pequena correção quadrática do regularizador. Para tornar o corte
renovável explícito, seja $nu_t^s$ o dual de
$g_t^s+c_t^s-overline(g)_t^s=0$. Em um ponto interior de geração,

$ -lambda_t+nu_t^s=0. $

Se o corte também é positivo,

$ pi_s+nu_t^s=0, $

e, sem o regularizador, $lambda_t=-pi_s$. Essa é a origem do piso de preço-sombra
em excesso renovável.

Se $mu_t$ multiplica
$S_(t+1)-S_t-A_t+h_t+v_t=0$, a estacionariedade de $h_t$ fornece

$ lambda_t=pi_h+mu_t $

em um ponto interior. A derivada em relação a $S_t$ transporta $mu_t$ entre horas;
quando o estoque toca piso ou teto, entram os multiplicadores das bandas. Assim,
$mu$ é o valor intertemporal da água armazenada, enquanto $pi_h$ é o custo de
oportunidade direto colocado no objetivo.

#keybox([O que o preço-sombra representa])[ 
  $lambda$ é o custo marginal do *modelo* para atender uma unidade adicional de
  carga no nó e na hora, dadas todas as restrições ativas. Ele incorpora mérito
  térmico, água, rampa, déficit e congestionamento. Não é o CMO divulgado, não é o
  PLD e não deve ser usado como se fosse uma observação de mercado.
]

== Impacto dos termos quadráticos e sensibilidade a $gamma_r$

Os três termos quadráticos têm papéis distintos. $epsilon_n n_b^2/2$ dá curvatura
a cada bloco térmico; $gamma_r(n_t-n_(t-1))^2$ pune variações horárias da térmica
agregada; e $epsilon_x||x||_2^2$ é Tikhonov sobre todo o vetor primal. A KKT acima
mostra exatamente como eles entram no preço marginal: $epsilon_n$ acrescenta
$epsilon_n n_b$, $gamma_r$ acrescenta uma segunda diferença temporal e Tikhonov
acrescenta um termo proporcional à variável. Eles não alteram apenas a velocidade
do solver.

#text(size: 8.15pt)[
#table(
  columns: (1.40fr, 1.05fr, 1.00fr, 1.00fr),
  inset: 3.7pt,
  align: right,
  [Cenário], [$lambda$ médio], [Rampa térmica RMS MW], [Vertim. TWh],
  [Baseline], [277,041], [40,44], [22,629],
  [Sem $epsilon_n$], [276,966], [47,48], [22,629],
  [Sem $gamma_r$], [277,031], [90,43], [22,629],
  [Sem Tikhonov], [272,950], [0,00007], [22,629],
  [Sem quadráticos], [272,908], [7,77], [22,629],
)
]

#imgfig(
  [Impacto dos termos quadráticos no despacho.],
  "outputs/figures/rel_quadratic_terms_sensitivity.png",
  note: [
    Os termos quadráticos são relevantes, mas para uma dimensão específica do
    resultado. A energia térmica permanece em 97,543 TWh e o vertimento em 22,629 TWh
    em todos os cenários na precisão reportada. Portanto, eles quase não alteram
    *quanto* cada recurso produz no ano. O que muda é *quando* a térmica produz e qual
    solução o QP seleciona entre despachos de custo linear muito parecido. Sem
    $gamma_r$, a rampa RMS mais que dobra, de 40,44 para 90,43 MW. Sem
    $epsilon_n$, ela sobe de forma mais moderada, para 47,48 MW. Sem Tikhonov, a
    solução selecionada distribui a térmica quase uniformemente e a rampa cai para
    perto de zero, ao mesmo tempo que o preço-sombra médio diminui cerca de
    4,09 R\$/MWh.

    Assim, os termos quadráticos não são relevantes para explicar os totais anuais
    deste caso, mas são muito relevantes para a forma horária, a unicidade prática da
    solução e o dual reportado. Tikhonov merece cautela especial: seu efeito mostra
    que havia muitos despachos quase equivalentes no objetivo linear e que o
    regularizador escolhe um deles; esse despacho suave não deve ser confundido com
    uma propriedade física estimada das usinas.
    O experimento está no notebook 07, §4.5.1, célula `[13]`; procure
    `QUADRATIC_SCENARIOS`, `quadratic_rows` e `quadratic_sensitivity`. A tabela é
    `dispatch_quadratic_terms_sensitivity.csv`; a figura está no notebook 08, §15,
    célula `[16]`.
  ]
)

Isolando $gamma_r$, a progressão $0$, $0,0005$, $0,001$, $0,002$, $0,005$ e $0,01$
reduz a rampa RMS de 90,43 para 19,21 MW — queda de 78,8% — enquanto o preço médio
varia apenas 0,0104 R\$/MWh e energia/vertimento permanecem iguais na precisão
reportada.

#imgfig(
  [Sensibilidade do despacho à penalidade de rampa $gamma_r$.],
  "outputs/figures/rel_gamma_r_sensitivity.png",
  note: [
    A sensibilidade a $gamma_r$ é grande para a variável que ele foi construído para
    controlar e pequena para os agregados do sistema. Quando $gamma_r$ passa de zero
    para 0,002, a rampa RMS cai de 90,43 para 40,44 MW, uma redução de 55%; em 0,01,
    chega a 19,21 MW. Em contraste, o preço-sombra médio varia apenas
    0,0104 R\$/MWh em toda a grade, e energia térmica, déficit e vertimento não mudam
    na precisão reportada. A curva azul, portanto, revela uma sensibilidade temporal
    forte; a laranja, uma sensibilidade econômica média quase nula.

    O valor $gamma_r=0,002$ produz suavização material, mas continua sendo uma escolha
    de regularização, não uma estimativa de custo físico de rampa. A figura permite
    avaliar o efeito numérico dessa escolha; não permite afirmar que 0,002 seja o
    parâmetro empiricamente correto. Resultados completos:
    `gamma_r_sensitivity.csv`, notebook 07, §4.5.1, célula `[13]`; procure
    `gamma_r_rows` e o laço `for gamma_r in [...]`. Reprodução gráfica: notebook 08, §15, célula
    `[16]`.
  ]
)

== Energia agregada não garante adequação horária

Mesmo que

$ sum_t (overline(g)_t^s+overline(g)_t^w+overline(n)_t+A_t)
  >= sum_t L_t, $

não segue que

$ overline(g)_t^s+overline(g)_t^w+overline(n)_t
  +overline(h)_t+F_t^"net,max" >= L_t quad "para todo" t. $

A primeira desigualdade compara energia acumulada. A segunda é uma condição
horária de potência e localização; ainda assim, é apenas necessária quando rampas,
estoques e fluxos acoplam horas e nós. Essa diferença formal explica por que um dia
com energia renovável abundante pode apresentar déficit no entardecer.

= Hidrologia, valor da água e sensibilidade

== Variação do custo de oportunidade da água

O clearing depende de um parâmetro que não observamos diretamente: quanto vale a
água guardada. Turbinar hoje é abrir mão de gerar amanhã, e esse custo de
oportunidade — o preço $pi_h$ que colocamos no objetivo por MWh hídrico — não sai
dos dados, é uma escolha. Como o resto do relatório usará esse valor da água para
precificar a solar, precisamos saber o quanto os resultados dependem dele antes de
confiar em qualquer número. É isso que este exercício faz: mantém a afluência de
referência e varia $pi_h$ em 80, 250 e 500 R\$/MWh, reportando a média de $lambda$
ponderada pelos pesos dos medoides e os totais anuais.

#text(size: 8.7pt)[
#table(
  columns: (1.02fr, 1.05fr, 1.12fr, 1.12fr, 1.12fr),
  inset: 4pt,
  align: right,
  [$pi_h$ R\$/MWh], [$lambda$ médio], [Hidro TWh], [Térmica TWh], [Vertimento TWh],
  [80], [241,53], [390,27], [80,23], [5,32],
  [250], [277,04], [372,96], [97,54], [22,63],
  [500], [504,27], [328,53], [141,98], [67,06],
)
]

Quando turbinar fica mais caro no objetivo, a hidro diminui, a térmica aumenta e o
preço-sombra sobe. Com fechamento cíclico e afluência fixa, parte da água não
turbinada precisa ser vertida; por isso vertimento também cresce. O exercício é
uma sensibilidade paramétrica, não uma calibração do verdadeiro valor da água do
SIN.

== Cenários seco, base e úmido

Para isolar hidrologia, multiplicamos todas as afluências equivalentes por 0,70,
1,00 e 1,30, mantendo demanda, rede, custos, capacidades e perfis renováveis
fixos.

#text(size: 8.25pt)[
#table(
  columns: (1.18fr, 0.92fr, 0.90fr, 0.90fr, 0.90fr, 0.90fr),
  inset: 3.6pt,
  align: right,
  [Cenário], [$lambda$ médio], [Térmica TWh], [Hidro TWh], [Vertim. TWh], [Déficit TWh],
  [Seco -30%], [7.474,26], [141,56], [276,02], [0,89], [52,93],
  [Base], [277,04], [97,54], [372,96], [22,63], [0,00],
  [Úmido +30%], [254,51], [93,94], [376,57], [137,70], [0,00],
)
]

O parâmetro $pi_h$ mede quanto *cobramos* pela água; falta ver o que acontece quando
muda a quantidade de água que de fato chega. Queremos, com o próximo exercício,
separar essas duas coisas e enxergar a não linearidade do sistema: um sistema
hidrotérmico não responde de forma suave à seca, ele tem um ponto de ruptura, e é
esse comportamento que a figura procura expor multiplicando todas as afluências por
0,70, 1,00 e 1,30 com o resto fixo.

#imgfig(
  [Sensibilidade hidrológica com as demais entradas fixas.],
  "outputs/figures/rel_hydrology_sensitivity.png",
  note: [
    Os três painéis mostram, da esquerda para a direita, o preço-sombra médio anual
    ponderado, o déficit anualizado e a decomposição entre térmica, hidro e
    vertimento. A leitura central é a *assimetria brutal* entre os cenários: tirar
    30% de afluência não encarece o sistema em 30%: o preço médio salta de 277,04
    para 7.474,26 R\$/MWh e surgem 52,93 TWh de déficit — enquanto adicionar 30%
    reduz o preço apenas para 254,51. Essa é a marca de um
    sistema que opera perto do limite de adequação: há muito a perder na seca e pouco
    a ganhar na cheia. Vale uma leitura fina de dois números que parecem contradizer
    a intuição. Primeiro, os 7.474,26 R\$/MWh ainda ficam *abaixo* do custo de
    déficit de 8.327 — não porque falte escassez, mas porque é uma média sobre nós e
    horas, e muitas horas seguem baratas; a média esconde os picos onde o preço bate
    o teto. Segundo, no cenário úmido o vertimento dispara para 137,70 TWh: água em
    excesso não vira valor, vira desperdício, porque o reservatório equivalente e o
    fechamento cíclico não deixam guardá-la para depois. As ressalvas pesam: o choque
    é uniforme e determinístico, não representa uma distribuição de energia natural
    afluente nem a dependência entre bacias, e uma seca real, concentrada em algumas
    regiões e persistente por meses, teria efeito distinto do que este corte parelho
    de 30% sugere.
  ]
)

== Por que o vertimento fica alto?

Com estoque cíclico em cada medoide, somar a equação do reservatório nas 24 horas
elimina a variação de estoque e impõe a identidade

$ sum_t A_(d,ell,t)=sum_t h_(d,ell,t)+sum_t v_(d,ell,t). $

Logo, toda afluência equivalente que não é turbinada precisa ser vertida no próprio
dia. A auditoria regional do baseline localiza 15,56 TWh no Norte (24,1% de sua
afluência equivalente) e 7,07 TWh no Sudeste/CO (3,1%); Sul e Nordeste são
numericamente nulos.

#text(size: 8.35pt)[
#table(
  columns: (0.75fr, 1.20fr, 1.15fr, 1.05fr),
  inset: 4pt,
  align: right,
  [Subs.], [Afluência TWh], [Vertim. TWh], [Vertim./afluência],
  [SE], [231,54], [7,07], [3,1%],
  [S], [73,39], [0,00], [0,0%],
  [NE], [26,05], [0,00], [0,0%],
  [N], [64,60], [15,56], [24,1%],
)
]

#imgfig(
  [Diagnóstico regional do vertimento no baseline.],
  "outputs/figures/rel_spill_diagnostic.png",
  note: [
    O valor alto é numericamente consistente, mas não deve ser interpretado como
    estimativa de vertimento físico do SIN. A entrada $A$ é uma proxy de energia
    hídrica observada distribuída no dia, não ENA física; $pi_h$ cobra apenas a água
    turbinada, enquanto verter é gratuito; e o fechamento diário impede transferir
    água entre medoides ou estações. Por isso o modelo pode preferir térmica marginal
    abaixo de $pi_h$ e verter a proxy restante. Dois testes confirmam o diagnóstico:
    baixar $pi_h$ de 250 para 80 R\$/MWh reduz o vertimento para 5,32 TWh; elevar a
    afluência em 30%
    leva-o a 137,70 TWh. O cálculo está no notebook 07, §4.5.1, célula `[13]`;
    procure `spill_by_loc` e a gravação de `spill_diagnostic_by_region.csv`. A
    figura está no notebook 08, §15, célula `[16]`.
  ]
)

#resultbox([Interpretação econômica da água])[ 
  O valor sistêmico de um MWh solar depende do estado hidrológico. Em escassez, a
  solar pode preservar água para horas sem sol e reduzir déficit; em abundância, o
  benefício marginal cai e o risco de vertimento cresce. Portanto, uma única
  receita média por MW não é estruturalmente invariável ao regime de água.
]

#pagebreak(weak: true)

= Laboratório situacional finito

O jogo de campo médio é o destino, mas ele é abstrato, e convém primeiro ver o
mecanismo de canibalização acontecer com atores concretos e contáveis. Por isso,
antes de tomar o limite de muitos agentes, montamos um laboratório com apenas dois
projetos solares no Sudeste/CO e perguntamos: o que acontece com o incumbente quando
um segundo projeto entra e cresce? O produtor 1 é o incumbente, com capacidade
$K_1 = 11.446,23$ MW. O produtor 2 é o entrante, cuja capacidade é uma fração $kappa$
da do incumbente, que varremos de zero a uma vez e meia,

$ K_2=kappa K_1, quad kappa in (0;0,25;dots;1,50). $

Ambos enfrentam exatamente o mesmo clearing da Seção 4 — com rede, térmica, hidro,
reservatório, eólica, nuclear e déficit —, de modo que qualquer efeito vem da
entrada, não de uma mudança de regras. A geração *potencial* de cada produtor $i$ é
$overline(g)_(i,d,t)^s = K_i xi_i a_(d,"SE",t)^s$, onde $xi_i$ é o fator de
produtividade do tipo do projeto e $a_(d,"SE",t)^s$ é o perfil solar do Sudeste/CO;
quando o nó precisa cortar, o corte é rateado proporcionalmente entre os dois. Com a
geração efetivamente entregue $g_(i,d,t)^s$ e o preço $lambda_(d,"SE",t)$ do
clearing, a receita anual do produtor (em R\$) e a receita capturada por unidade de
energia entregue (em R\$/MWh) são

$ Rev_i=sum_d p_d sum_t lambda_(d,"SE",t) g_(i,d,t)^s, $

$ R_i^"cap"=frac(Rev_i, sum_d p_d sum_t g_(i,d,t)^s). $

O lucro desconta o custo anualizado da capacidade, $q_("SE",theta_i)$ em R\$/MW/ano
para o tipo $theta_i$ do projeto:

$ Pi_i=Rev_i-q_("SE",theta_i)K_i. $

Um detalhe de honestidade numérica: quando $kappa=0$ o entrante simplesmente não
existe, então sua receita é zero e sua receita capturada por MWh é *indefinida* — o
denominador também é zero. Preferimos deixar a lacuna a preenchê-la com um zero
artificial que distorceria a leitura da curva.

Queremos ver, na figura, duas coisas ao mesmo tempo: se a entrada corrói a receita
do incumbente (a canibalização) e se essa corrosão é gradual ou tem saltos. A
resposta ao segundo ponto é a mais reveladora e não é óbvia à primeira vista.

#imgfig(
  [Incumbente e entrante: receita capturada, lucro e curtailment.],
  "outputs/figures/rel_finite_two_player.png",
  note: [
    No eixo horizontal está $kappa=K_2/K_1$; os painéis mostram, em sequência, a
    receita capturada por MWh, o lucro anual e o corte alocado a cada projeto. No
    primeiro painel, a receita capturada é a receita total dividida pela energia
    efetivamente entregue por cada agente. As duas curvas coincidem quando ambos
    produzem com o mesmo perfil e recebem o mesmo rateio do corte. À
    medida que o entrante cresce, essa receita por MWh cai em degraus, de cerca de
    319 para 255 R\$/MWh, mostrando que os MW adicionais reduzem o valor médio da
    produção simultânea.

    O painel central torna a diferença entre os papéis mais concreta. O incumbente
    começa com lucro anual alto e perde espaço conforme a concorrência aumenta; o
    entrante parte de zero, ganha escala e atinge seu melhor resultado perto de
    $kappa=1$, mas não cresce indefinidamente. Em $kappa=1,25$, seu lucro recua, apesar
    da capacidade maior, porque a queda da receita capturada e o custo da capacidade
    passam a pesar mais que a energia adicional. Em $kappa=1,50$, há recuperação
    parcial, sinal de que o clearing por blocos produz quinas e não uma resposta
    perfeitamente lisa.

    O corte anualizado no painel direito é da ordem de $10^(-6)$ GWh, isto é,
    desprezível. A queda de rentabilidade não é causada por energia solar impedida de
    entrar, mas pela redução do valor capturado e pelo custo de instalar capacidade.
    O gráfico mostra, de forma simples, que a entrada beneficia o entrante até certo
    ponto e, ao mesmo tempo, corrói a posição do incumbente; nenhum dos dois efeitos
    cresce de forma linear porque cada escala do entrante resolve um despacho com
    restrições ativas potencialmente diferentes.
  ]
)

#keybox([Por que este laboratório é uma ponte e não o MFG])[ 
  Os dois produtores alteram deliberadamente a capacidade total em uma análise de
  sensibilidade. No MFG, cada agente individual é infinitesimal e toma o campo como
  dado; apenas a distribuição agregada altera o clearing. O laboratório revela o
  mecanismo de canibalização e de mudança de restrições ativas, mas não substitui
  a condição de consistência populacional.
]

= Investimento estático e planejador

== Custos anualizados e tipos

O custo base por MW combina CAPEX, conexão e O&M por meio do fator de recuperação
de capital. O tipo `high` tem produtividade relativa 1,08 e multiplicador de custo
0,94; o tipo `low` tem produtividade 0,94 e multiplicador de custo 1,08. As massas
iniciais são 45% e 55%, respectivamente.

Precisamos agora de um objeto que ligue *onde e como o parque está distribuído* à
receita que cada tipo de projeto obtém. Chamamos esse objeto de $R_(ell,theta)(m)$: a
receita spot anual por MW de um projeto na localização $ell$ e tipo $theta$, dada uma
distribuição $m$ do parque. Aqui $m$ é a peça nova e central — é o *campo*, isto é, a
descrição de como a capacidade está espalhada entre localizações e tipos; é o objeto
que o jogo de campo médio vai tornar endógeno. A receita é

$ R_(ell,theta)(m)=sum_d p_d sum_t
  lambda_(d,ell,t)(m) a_(d,ell,t)^s xi_theta
  delta_(d,ell,t)^"delivery"(m), $

somando, sobre medoides (peso $p_d$) e horas, o preço $lambda(m)$ que o clearing
devolve, multiplicado pela energia que um MW daquele tipo gera — o perfil
$a^s$ vezes a produtividade $xi_theta$ — e pela fração de fato entregue, que
desconta o curtailment. Para fechar também as horas sem recurso solar, o notebook
implementa explicitamente a convenção

$
  delta_(d,ell,t)^"delivery"=cases(
    frac(g_(d,ell,t)^s, overline(g)_(d,ell,t)^s)
      & "se " overline(g)_(d,ell,t)^s > 10^(-6),
    1 & "se " overline(g)_(d,ell,t)^s <= 10^(-6),
  ), quad 0 <= delta_(d,ell,t)^"delivery" <= 1.
$

#nbref([
  Notebook 07, §4.2, célula `[11]`, método `DispatchEngineCVX.solve`: procure
  `delivery = np.clip(...)` e, logo depois, o cálculo de `R[(s, agent_type)]`,
  que multiplica `lambda`, `cf_agent` e `delivery`. Essa é a implementação direta
  da equação acima.
])

O segundo ramo evita a indeterminação $0/0$ à noite. Escolher 1 nesse ramo — em
vez de 0 — é a convenção numérica efetivamente usada (`out=ones_like` antes do
truncamento ao intervalo $[0,1]$); ela não cria receita artificial, pois nessas horas
$a_(d,ell,t)^s=0$ e, portanto, o produto que entra em $R$ continua nulo.
A dependência em $m$ não é decorativa, é o cerne do problema: a distribuição $m$
determina a capacidade instalada, a capacidade alimenta o clearing, e o clearing fixa
o preço e, portanto, a própria receita $R(m)$. É esse laço — a receita depende da
distribuição que a receita ajuda a formar — que impede tratar o investimento como uma
decisão isolada e exige um ponto fixo.

A canibalização foi derivada na teoria e medida no CMO observado; falta confirmá-la
no operador completo, com rede, água e déficit todos ativos. A figura seguinte serve
a isso — traça $R$ contra capacidade — e, de quebra, mostra se a relação é bem
comportada o bastante para que o ponto fixo de investimento tenha chance de convergir.

#imgfig(
  [Receita anual por MW como função da capacidade solar agregada.],
  "outputs/figures/canibalizacao_R_de_C.png",
  note: [
    O eixo horizontal escala ao mesmo tempo a capacidade solar de referência das
    quatro regiões e cada curva é a receita anual por MW de um projeto em uma
    localização, recalculada do zero pelo clearing a cada ponto. Todas as curvas
    descem, o que fecha o argumento: a canibalização derivada com o CMO observado
    reaparece no operador físico completo, não era artefato da regressão. Mas as
    curvas não são paralelas, e a diferença entre elas é a mensagem locacional —
    níveis e inclinações separam-se por rede, perfil renovável e recurso disponível,
    de modo que o mesmo MW adicional vale coisas distintas conforme onde entra, e o
    ritmo com que a receita se esgota também difere. Uma cautela metodológica sincera:
    a curva cair de forma monótona nessa grade é uma boa notícia para a estabilidade
    do ponto fixo — receita que cai com entrada empurra o sistema de volta ao
    equilíbrio —, mas monotonicidade observada numa grade não *prova* que o operador
    de melhor resposta seja uma contração global, então a convergência que veremos a
    seguir é verificada numericamente, não garantida de antemão.
  ]
)

== Benchmark estático de campo médio

Com a receita $R(m)$ em mãos, o passo natural é perguntar quanta capacidade um
investidor racional escolheria. O benchmark estático responde a isso na versão mais
simples possível: a massa de projetos em cada região fica *fixa* (ninguém migra) e a
única decisão é *quanto* investir. Cada investidor de tipo $theta$ na localização
$ell$ escolhe um índice de capacidade $k$ para maximizar seu ganho líquido

$ U_(ell,theta)(k)
  = (R_(ell,theta)-q_(ell,theta)+p) kappa_u k
  - frac(gamma,2) k^2. $

O primeiro termo é o lucro por MW — a receita $R_(ell,theta)$ menos o custo
anualizado $q_(ell,theta)$, mais um prêmio não spot opcional $p$ que representa
receita fora do mercado de curto prazo (contratos) — convertido em MW pelo fator
$kappa_u$ e multiplicado pela capacidade $k$. Nesta execução,
$kappa_u="39.069,292"/160="244,183"$ MW por bin: $k$ é um índice adimensional e
$kappa_u k$ é sua capacidade física. Esse valor é apenas a normalização que faz o
bin inicial agregado reproduzir a proxy solar; não é o degrau de 93 MW usado mais
adiante no MFG dinâmico. Mais precisamente, 160 é `K0_BIN`: o índice de referência
em que toda a massa é colocada na inicialização do benchmark estático. Como a massa
normalizada soma um, $kappa_u times 160="39.069,292"$ MW recupera exatamente a
capacidade agregada observada no ponto de partida.

A grade completa usada no notebook é

$ k in cal(K)^"static"={0,1,dots,480}, quad abs(cal(K)^"static")=481. $

#nbref([
  Notebook 07, “Parâmetros de velocidade e do modelo”, célula `[2]`: procure
  `K0_BIN = 160`, `K_MAX_BIN = 3 * K0_BIN`, `K_BINS = np.arange(...)` e
  `N_K`. A conversão do índice em MW, `CAP_UNIT_MW`, é calculada no notebook 07,
  §2, célula `[7]`; a melhor resposta que usa essa grade está em §6.1, célula
  `[16]`, função `best_response`.
])

Assim, o limite superior corresponde a três vezes o índice de referência
($480=3 times 160$), ou 117,208 GW agregados caso toda a massa estivesse no último
bin. O segundo termo, $frac(gamma,2)k^2$, é
uma penalidade quadrática crescente que representa custos marginais de expansão que
sobem com o tamanho: $gamma$ controla essa curvatura e é, na prática, o que impede
que o investidor queira capacidade infinita quando o lucro por MW é positivo.
Maximizando, a capacidade ótima interior é

$ k^* = pos(frac((R-q+p)kappa_u,gamma)), $

onde $pos(dot)$ zera valores negativos (não se investe quando o lucro por MW é
negativo), truncada explicitamente ao intervalo $[0,480]$. Quando $k^*$ cai entre dois nós $b$ e
$b+1$ da grade discreta, ele é representado por uma média ponderada dos dois, com
peso $phi=k^*-b$ no nó superior; essa representação baricêntrica preserva
exatamente a capacidade média, e é importante frisar que ela *não* é uma "mistura
ótima" de um jogo discreto — os dois bins vizinhos não precisam dar o mesmo ganho.

O equilíbrio é um *ponto fixo*: a distribuição de capacidade que os investidores
escolhem tem de ser a mesma que gerou as receitas às quais eles responderam. Escrito
com o operador de melhor resposta $BR$ (que, dado o campo de receitas, devolve a
distribuição $k^*$ ótima),

$ m=BR(R(m)). $

Em palavras: parte-se de uma distribuição $m$, calcula-se a receita $R(m)$ pelo
clearing, cada agente responde com sua capacidade ótima formando $BR(R(m))$, e
procura-se o $m$ em que essa resposta reproduz o próprio $m$.

A rotina numérica usa relaxação constante $alpha=0,06$:

$ m^(n+1)=(1-alpha)m^n+alpha BR(R(m^n)). $

O tamanho do passo

$ s_m^n=||m^(n+1)-m^n||_infinity $

não certifica equilíbrio, pois pode ser pequeno apenas porque $alpha$ é pequeno.
Os resíduos válidos são calculados antes da relaxação:

$ r_m=||BR(R(m))-m||_infinity, $

$ r_C=max_ell frac(abs(C_ell(BR)-C_ell(m)),C_ell(m)+1), $

$ r_R=max_(ell,theta)
  frac(abs(R_(ell,theta)(BR)-R_(ell,theta)(m)),abs(R_(ell,theta)(m))+1). $

O benchmark só é rotulado como certificado quando os três resíduos ficam abaixo
das respectivas tolerâncias. Nesta execução, os limites são

$ r_m < 5 times 10^(-4), quad r_C < 5 times 10^(-3), quad r_R < 5 times 10^(-3). $

Essas desigualdades reproduzem o teste estrito implementado no notebook: não basta
que a iteração relaxada dê passos pequenos; os três resíduos do operador sem
relaxação precisam passar simultaneamente.

#nbref([
  Notebook 07, §6.2 “Ponto fixo e certificado do benchmark estático”, célula
  `[17]`, função `solve_mfg`. No corpo da iteração, procure
  `residual_m = ...` e `residual_C = ...`; depois do clearing final, procure
  `residual_R = ...`. A decisão conjunta está no bloco
  `certified = bool(residual_C < TOL_C_FP and residual_m < TOL_M_FP and
  residual_R < TOL_R_FP)`. A execução e a impressão dos três valores ficam em
  §6.4, célula `[19]`.
])

== Calibração por entrada livre

O nível agregado é calibrado para a proxy observada de 39,069 GW. Com a rede
a rotina seleciona como objeto de fechamento o prêmio não spot: mantém
$gamma=0,312$ milhão de R\$ por bin ao quadrado e calibra
$p=19,1$ mil R\$/MW/ano. A solução final contém 38,954 GW, diferença de 0,115 GW
em relação ao alvo dentro das tolerâncias do ponto fixo.

#nbref([
  Notebook 07, §6.3, célula `[18]`: `calibrate_gamma`,
  `calibrate_premium` e `calibrate_entry_wedge`. A chamada que escolhe o objeto,
  executa `solve_mfg` e grava `entry_wedge_calibration_path.csv` está em §6.4,
  célula `[19]`.
])

#warnbox([O que a calibração não prova])[ 
  Reproduzir o total observado é uma condição de calibração, não uma validação fora
  da amostra. A distribuição regional, a sensibilidade a custos, os preços e a
  dinâmica continuam sendo resultados condicionais. Como a massa regional fica
  fixa nessa camada, ela escolhe capacidade dentro da localização atribuída; a
  migração endógena pertence ao MFG dinâmico.
]

Afirmar que existe equilíbrio exige mostrar que a iteração de fato converge e — mais
sutil — que ela convergiu para um ponto fixo *verdadeiro*, não apenas parou de se
mexer porque o passo era pequeno. A figura serve para separar essas duas coisas, uma
distinção metodológica que o relatório leva a sério.

#imgfig(
  [Convergência e capacidade do benchmark estático.],
  "outputs/figures/rel_mfg_convergencia.png",
  note: [
    À esquerda, em escala logarítmica, estão os resíduos do operador ($r_C$, $r_m$) e
    o tamanho do passo relaxado; à direita, a capacidade agregada caminhando para o
    alvo de calibração. A rotina encerra na iteração 132 com
    $r_C=3,37 times 10^(-4)$, $r_m=4,88 times 10^(-4)$ e $r_R=4,00 times 10^(-7)$, os
    três abaixo das tolerâncias, e por isso o ponto é rotulado como certificado. O
    detalhe que merece atenção — e que é uma armadilha comum neste tipo de trabalho —
    é que a curva do passo fica *sistematicamente abaixo* dos resíduos, por um fator
    ligado ao amortecimento $alpha$. Se tivéssemos usado o passo como critério de
    parada, declararíamos convergência cedo demais e com folga falsa; é justamente
    por isso que a certificação se apoia nos resíduos do operador *sem* relaxação, e
    não no passo. A convergência aqui é uma propriedade verificada do algoritmo, não
    uma prova de que o equilíbrio seja único.
  ]
)

#text(size: 8.7pt)[
#table(
  columns: (0.72fr, 1.00fr, 1.00fr, 1.10fr),
  inset: 4pt,
  align: right,
  [Subs.], [Proxy observada GW], [Privado GW], [Privado - proxy GW],
  [SE], [19,077], [18,749], [-0,328],
  [S], [5,919], [2,698], [-3,221],
  [NE], [11,789], [16,234], [+4,445],
  [N], [2,284], [1,273], [-1,011],
  [*Total*], [*39,069*], [*38,954*], [*-0,115*],
)
]

A calibração aproxima o total, mas o campo de receitas redistribui capacidade dentro
das massas regionais: o Nordeste cresce e as demais regiões ficam abaixo das
proxies. Isso é uma implicação interna do modelo com a ligação NE--SE, não uma
recomendação de localização.

== Contrafactual de total fixo

Para preservar a ideia de “orçamento fixo” sem chamá-la de equilíbrio privado,
resolvemos o problema do planejador com a restrição

$ sum_ell K_ell=K^"obs"=39,069 " GW". $

Capacidade total e todas as demais entradas ficam fixas; apenas a localização é
otimizada junto com o despacho. O resultado distribui 19,327 GW no Sudeste/CO e
19,742 GW no Nordeste, com valores numericamente nulos no Sul e no Norte. Essa
solução de canto mostra que, sob os custos e a rede reduzida, o operador valoriza
o corredor SE--NE. Ela deve ser lida como diagnóstico de identificação: uma
alocação tão concentrada sinaliza que custos locacionais, limites de rede e
potenciais máximos ainda precisam de calibração externa antes de uso normativo.

== Planejador social com total livre e limites regionais

O investidor privado enxerga só a própria receita; o planejador social enxerga o
sistema inteiro. A diferença entre os dois é o coração da pergunta normativa — o
mercado investe de menos ou de mais do ponto de vista social? — e por isso resolvemos
também o problema do planejador, que escolhe capacidade $K$ e despacho $x$ de uma vez
para minimizar o custo total do sistema:

$
  min_(K,x) sum_ell q_ell K_ell + C^"op"(K,x),
$

onde o primeiro termo é o custo anualizado de instalar a capacidade ($q_ell$ por MW)
e $C^"op"(K,x)$ é o custo operacional que o próprio clearing gera — combustível
térmico, água, déficit e cortes. Tudo sujeito ao mesmo clearing da Seção 4 e a um
teto regional que impede soluções extremas,

$ 0<=K_ell<=3K_ell^"obs". $

#nbref([
  Notebook 07, §8, célula `[22]`, função `solve_social_planner`. Procure
  `K = cp.Variable`, o termo `objective += q_base[s] * K[i]`, a restrição
  `K[i] <= 3.0 * K_SOLAR_OBS[s]` e a chamada
  `cp.Problem(cp.Minimize(objective), constraints)`.
  A mesma célula resolve o caso livre (`social`) e o total fixo (`fixed_budget`) e
  grava `private_vs_social_capacity.csv`. A figura correspondente está no notebook
  08, §9, célula `[10]`.
])

Para tornar a comparação auditável, é preciso registrar que os dois problemas não
diferem apenas por “internalizar uma externalidade”. Eles compartilham o mesmo
núcleo de clearing operacional, mas não a mesma agregação de tipos, a mesma função de investimento nem as mesmas
restrições de escolha:

#text(size: 7.7pt)[
#table(
  columns: (1.12fr, 1.75fr, 1.75fr, 1.55fr),
  inset: 3.2pt,
  align: left,
  [Elemento], [Privado usado na figura], [Planejador social], [Implicação],
  [Objetivo], [Maximiza receita própria menos custo por tipo e curvatura $gamma k^2/2$.], [Minimiza custo anualizado de capacidade mais custo operacional sistêmico.], [Receita spot é ganho privado, mas transferência no problema social.],
  [Custo linear], [$q_(ell,theta)=q_ell^"base" times (0,94 " ou " 1,08)$.], [$q_ell^"social"=q_ell^"base"$, sem dimensão de tipo.], [O custo social de referência é cerca de 1,7% menor que a média privada ponderada.],
  [Produtividade solar], [$a_(d,ell,t)^s sum_theta xi_theta K_(ell,theta)$, com $xi_theta$ por tipo.], [$a_(d,ell,t)^s K_ell$, produtividade unitária agregada.], [O planejador também colapsa a heterogeneidade produtiva.],
  [Curvatura estática $gamma$], [Entra e é calibrada para reproduzir 39,069 GW.], [Não entra.], [Parte do hiato decorre dessa diferença de especificação.],
  [Entrada, mudança e concentração], [Não entram neste benchmark estático; entram apenas no MFG dinâmico.], [Não entram.], [Não explicam diretamente as barras privada e social da figura.],
  [Capacidade e localização], [Grade de bins e massas regionais fixas.], [$K_ell$ contínuo e livre, sujeito a $K_ell<=3K_ell^"obs"$.], [O planejador pode realocar e expandir com mais liberdade.],
  [Operação], [Receita por MW calculada pelo clearing e tomada como dada pelo agente infinitesimal.], [Mesmo custo operacional, escolhido conjuntamente com $K$.], [O planejador internaliza térmica, água, déficit, corte e rede.],
)
]

O valor $q_ell^"base"$ não vem de um retorno social estimado: ele é construído como
$("CAPEX"_ell+"conexão"_ell) "FRC"+"O&M"_ell$, com WACC de 8% e vida de 25 anos.
Como o planejador escolhe apenas $K_ell$ agregado, sem $theta$, adotamos esse custo
de referência como $q_ell^"social"$. Já a média privada dos multiplicadores vale
$0,45 times 0,94+0,55 times 1,08=1,017$. Portanto, até o custo linear não é
exatamente igual entre os problemas. A escolha é transparente e computacionalmente
coerente com a agregação social, mas impede atribuir todo o hiato de capacidade a
uma externalidade econômica.

A curvatura $gamma$ também foi omitida deliberadamente do planejador porque, no
benchmark privado, ela é uma forma reduzida calibrada para limitar a expansão
individual e fechar o total observado; não foi medida como custo real de terreno,
conexão ou cadeia de suprimento para a sociedade. Transportá-la mecanicamente para o
objetivo social daria interpretação de recurso a um calibrador privado. A decisão
evita essa falsa equivalência, mas deixa o planejador mais propenso a expandir e,
portanto, amplia potencialmente o hiato. De modo análogo, entrada, mudança e
concentração só aparecem na recompensa do MFG dinâmico. Custos fixos e de mudança
poderiam consumir recursos sociais na realidade, enquanto a concentração pode
representar tanto externalidade física quanto regularização; a execução corrente
não os decompõe nem os inclui no planejador. Essa assimetria deve ser harmonizada em
um exercício próprio antes de uma conclusão de bem-estar.

A diferença essencial em relação ao privado aparece no benefício marginal de mais um
MW, que pelo envelope é

$ frac(partial C^"op",partial K_ell)
  =-sum_d p_d sum_t lambda_(d,ell,t) a_(d,ell,t)^s
    delta_(d,ell,t)^"delivery" + "termos de corte". $

O planejador iguala o custo marginal de investir à redução marginal do custo
sistêmico. E é aqui que os dois divergem de forma conceitualmente limpa: o privado
toma $R(m)$ como dado e para quando *sua* receita cobre *seu* custo; o planejador
internaliza tudo o que um MW a mais faz ao sistema — menos déficit, menos
congestionamento, água preservada, menos térmica cara —, benefícios que o preço não
devolve inteiramente ao investidor individual.

A figura a seguir é uma das mais importantes do relatório: ela põe lado a lado as
quatro respostas para "quanta solar e onde" — a observada, a de orçamento fixo, a
privada e a social — para que a divergência entre mercado e planejador fique
visível de uma vez. É aqui que queremos ler se o modelo aponta subinvestimento e onde.

#imgfig(
  [Capacidade observada, total fixo, privada e social por subsistema.],
  "outputs/figures/rel_mercado_vs_social.png",
  note: [
    As barras cinza, laranja e verde são a proxy `p99`, o contrafactual de total
    fixo e o benchmark privado calibrado, todas próximas de 39,1 GW por construção; a
    azul é o planejador social, que escolhe 80,76 GW — pouco mais que o *dobro*. O
    resultado é grande demais para ser lido ingenuamente, e três leituras críticas
    são obrigatórias. Primeiro, as três primeiras barras são ancoradas no mesmo total,
    mas o salto do planejador combina benefícios sistêmicos internalizados com
    diferenças de especificação já documentadas — sobretudo ausência da curvatura
    privada, custo sem tipo e maior liberdade locacional. Sem um contrafactual com
    funções e restrições harmonizadas, não é possível decompor quantos GW do hiato
    vêm de externalidade e quantos vêm dessas escolhas. Segundo, nenhum teto regional
    $3K^"obs"$ está ativo nesta execução: SE, S, NE e N ficam em 51,23, 5,57, 22,06
    e 1,90 GW. Os limites ainda definem o domínio do problema, mas não determinam
    diretamente esta solução. Terceiro, nesta parametrização, o benefício social
    de mais solar (sobretudo por adequação)
    supera seu custo bem além do que o privado investe. Mas "nesta parametrização"
    carrega muito: o resultado é sensível ao custo de déficit, ao custo anualizado,
    ao teto e à ausência de contratos. A leitura correta não é "devemos dobrar a
    capacidade", e sim "há um fosso entre valor privado e social grande o bastante
    para justificar investigar por quê".
  ]
)

O contrafactual de total fixo praticamente zera Sul e Norte porque o planejador
precisa manter 39,069 GW, mas pode escolher onde colocá-los. O Nordeste reúne o menor
custo anualizado de referência, 437,3 mil R\$/MW/ano, boa produtividade solar e uma
interface direta de 8 GW com o Sudeste/CO. O Sudeste/CO, por sua vez, concentra a
maior carga e ocupa a posição central da rede reduzida; instalar ali evita depender
de uma exportação adicional para atender o principal centro consumidor. O Sul tem
custo maior que SE e NE e é um nó terminal conectado apenas por SE--S. O Norte
combina o maior custo regional, 514,9 mil R\$/MW/ano, com uma interface periférica de
4,5 GW pelo hub. Dentro desses custos e limites, um MW deslocado para S ou N reduz
menos o custo operacional do que o mesmo MW em SE ou NE, e a solução convexa vai para
o canto.

Esse quase zero é um resultado do laboratório, não uma recomendação para retirar
capacidade dessas regiões. O planejador não recebe benefício por diversificação
geográfica, segurança diante de extremos, disponibilidade de terra, limites de
conexão locais ou expansão futura da transmissão; também não impõe piso regional.
Na realidade, esses elementos dão valor a uma carteira mais espalhada. Portanto, a
concentração no corredor SE--NE é útil justamente por revelar quais mecanismos ainda
precisam ser calibrados ou acrescentados antes de uma leitura normativa.

#text(size: 8.35pt)[
#table(
  columns: (0.70fr, 0.95fr, 1.10fr, 0.95fr, 1.05fr),
  inset: 3.8pt,
  align: right,
  [Subs.], [Proxy GW], [Total fixo GW], [Privado GW], [Social limitado GW],
  [SE], [19,077], [19,327], [18,749], [51,234],
  [S], [5,919], [0,000], [2,698], [5,568],
  [NE], [11,789], [19,742], [16,234], [22,058],
  [N], [2,284], [0,000], [1,273], [1,904],
  [*Total*], [*39,069*], [*39,069*], [*38,954*], [*80,764*],
)
]

#resultbox([Leitura correta da comparação])[ 
  A comparação não identifica empiricamente “subinvestimento” ou “sobreinvestimento”.
  Ela mostra que, dentro do modelo calibrado, o agente privado e o planejador
  valorizam margens diferentes e enfrentam especificações diferentes. Sem repetir o
  exercício com custo, curvatura e limites harmonizados, os 41,810 GW de diferença
  não podem ser chamados integralmente de externalidade. O resultado social é
  sensível ao VOLL $pi_u$ definido na Seção 4, ao custo
  anualizado, ao teto regional e à representação de contratos. O fato de um teto
  nenhum teto está ativo nesta execução, mas sua forma e amplitude continuam
  condicionando o conjunto factível e devem ser testadas antes de qualquer conclusão
  normativa.
]

= Mean Field Game dinâmico locacional

== Por que o benchmark estático não basta

O benchmark anterior escolhe capacidade com localização fixa. Um MFG dinâmico deve
representar simultaneamente:

- o estado individual de capacidade, localização e tipo;
- uma política ótima condicionada ao campo agregado;
- a evolução forward da distribuição populacional;
- o clearing que transforma essa distribuição em capacidade, despacho e receita;
- uma condição de consistência sobre a trajetória inteira.

#keybox([Mapa da notação usada a partir desta seção])[
  O índice $y=0,dots,H$ identifica o estágio da decisão dinâmica, enquanto $n$
  identificará a iteração numérica usada para encontrar o ponto fixo. O estado
  individual é $x=(ell,theta,b)$: localização $ell$, tipo solar $theta$ e bin de
  capacidade $b$. A letra $m_y$ é a distribuição candidata da população no estágio
  $y$; $hat(m)_y$ é a distribuição que a política ótima implica quando fazemos a
  população avançar pela equação forward.

  A letra $K_(ell,y)$ é a capacidade em MW que resulta de $m_y$, e $hat(K)_(ell,y)$
  é a capacidade em MW que resulta de $hat(m)_y$. A letra maiúscula
  $R_(y,ell,theta)$ é a receita anual por MW, em R\$/MW/ano, devolvida pelo clearing
  sob $m$; $hat(R)$ é a receita recalculada depois de substituir $m$ por $hat(m)$.
  A política $pi_y(a|x)$ é a probabilidade de o estado $x$ escolher a ação $a$.

  Finalmente, $r_m$, $r_C$ e $r_R$ são resíduos numéricos: medem, respectivamente,
  a inconsistência da distribuição, da capacidade e da receita. O subscrito $C$ de
  $r_C$ é uma nomenclatura histórica para “capacidade”, embora as equações deste
  relatório usem a letra $K$ para capacidade. Resíduo pequeno significa que a
  resposta dos agentes quase reproduz o campo que eles tomaram como dado; resíduo
  zero significaria reprodução exata.
]

O horizonte dinâmico tem seis *decisões* ($y=0,dots,5$) e, por isso, sete
distribuições populacionais ($m_0,dots,m_6$). Um estágio é um ciclo abstrato de
revisão de investimento, não uma hora do despacho e não foi identificado como um
ano civil. As receitas que entram em cada recompensa já estão anualizadas pelos
pesos dos medoides e são contabilizadas uma vez por estágio. Consequentemente,
$beta=0,95$ é desconto *por estágio* do laboratório. Interpretá-lo como taxa anual
de 5% só seria válido sob a hipótese adicional de que um estágio equivale a um ano,
hipótese que esta execução não impõe. Assim, $H=6$ serve para estudar uma sequência
curta de ajustes e realocações, não para datar uma previsão em seis anos.

Vale explicar por que isso não é apenas o benchmark anterior "com mais passos". Num
jogo finito com muitos projetos, o agente $i$ maximizaria seu ganho de longo prazo
$max_(a_i) J_i(a_i,a_(-i))$, onde $a_i$ é a sua própria estratégia e $a_(-i)$ é o
vetor com as estratégias de *todos os outros* agentes — um objeto gigantesco e
impossível de rastrear. A ideia de campo médio é substituir esse vetor pela
*distribuição* $m_y$ da população no estágio $y$: em vez de saber o que cada rival
faz, o agente só precisa saber como a população está espalhada. Cada agente é então
infinitesimal — sozinho não move nada — e toma o par campo-receita $(m,R)$ como dado
ao resolver seu problema de decisão. O fecho do argumento, e a dificuldade real, é a
consistência: a população agregada de agentes que respondem a $(m,R)$ tem de produzir
de volta exatamente o mesmo $(m,R)$ que eles anteciparam.

== Estado, ação e transição

Para operacionalizar isso precisamos dizer o que descreve um agente, o que ele pode
fazer e como isso muda seu estado. O *estado* de um agente é o trio

$ x=(ell,theta,b) in cal(X)=cal(L) times cal(Theta) times cal(B): $

sua localização $ell$, seu tipo solar $theta$ (aquele `high`/`low` permanente da
Seção 2) e seu nível de capacidade $b$, um índice discreto de 0 a 12 em que cada
degrau vale 93 MW. Aqui $b$ é um bin de capacidade, distinto do índice dos blocos
térmicos usado no clearing. A *ação* combina para onde ir e quanto crescer,

$ a=(ell',delta), quad delta in (-1,0,1), $

isto é, escolher a localização de destino $ell'$ e decidir se reduz um degrau de
capacidade, mantém, ou sobe um degrau ($delta=-1,0,+1$). A transição é
determinística — não há incerteza sobre o resultado da ação —,

$ T((ell,theta,b),(ell',delta))
  =(ell',theta,b'), quad b'=b+delta, $

válida apenas para ações que mantêm $b'$ entre 0 e 12; o tipo $theta$ nunca muda. A
população é descrita pela distribuição $m_y$, que é uma medida de probabilidade
legítima — não negativa e somando um sobre todos os estados,

$ m_y(ell,theta,b)>=0, quad sum_(ell,theta,b)m_y(ell,theta,b)=1. $

Assim, $m_y(ell,theta,b)$ é literalmente a fração da população que, no estágio $y$,
está na localização $ell$, possui o tipo $theta$ e ocupa o bin $b$. Não é MW, não é
número absoluto de usinas e não é fluxo de energia: é uma participação entre zero e
um. Somar $m_y$ sobre todos os estados produz um; somá-la apenas sobre $theta$ e $b$
produz a participação populacional de uma localização.

A distribuição inicial é fixada por uma única expressão:

$
  m_0(ell,theta,b)=s_ell^"obs" w_theta bold(1)_{b=4}, quad
  s_ell^"obs"=frac(K_ell^"obs", sum_j K_j^"obs"), quad
  w_"high"=0,45, quad w_"low"=0,55.
$

Aqui $K_ell^"obs"$ é a capacidade solar observada usada na proxy e
$bold(1)_{b=4}$ vale um somente no bin 4. Portanto, todos os agentes equivalentes
começam com 372 MW, as massas entre subsistemas reproduzem as participações da
capacidade observada, e dentro de cada subsistema 45% da massa é do tipo `high` e
55% do tipo `low`. Como $sum_ell s_ell^"obs"=1$ e
$w_"high"+w_"low"=1$, essa fórmula também deixa explícito que $m_0$ soma um e é
mantido fixo durante a solução do ponto fixo.

Como $m_y$ é uma fração, precisamos de uma escala para traduzi-la em MW. Chamando $M$
o número equivalente de agentes, a capacidade instalada em cada nó é

$ K_(ell,y)(m_y)=M sum_(theta,b) K_b m_y(ell,theta,b),
  quad K_b=93b " MW", $

ou seja, a soma da capacidade $K_b$ de cada bin ponderada pela fração de agentes ali.
O valor de $M$ é fixado uma única vez, para que o estado inicial — concentrado em
$b=4$, isto é, 372 MW por agente — reproduza os 39,069292 GW da proxy observada:

$ M=frac("39.069,292 MW", 4 times 93 " MW")="105,025". $

Portanto, $M$ é o número de *agentes equivalentes* usado para converter frações em
MW. Não afirma que existam 105 usinas ou empresas reais, e tampouco é uma massa que
possa crescer; entrada e saída são movimentos entre o bin zero e bins positivos
dentro dessa população normalizada. A identidade ancora a escala inicial, mas não
impõe nada sobre a trajetória posterior.

== Recompensa individual

O que faz o agente escolher uma ação em vez de outra é sua recompensa a cada estágio.
Ela precisa capturar não só o lucro, mas também os atritos reais de investir — mudar
de tamanho custa, entrar custa, mudar de lugar custa, e aglomerar-se onde todos já
estão custa. Dada a receita $R_y$ que o clearing devolve, a recompensa de ir de
$(ell,b)$ para $(ell',b')$ é

$
  r_y(x,a;m_y,R_y)
  = (R_(y,ell',theta)+p-q_(ell',theta))K_(b')
  - frac(chi,2)(K_(b')-K_b)^2
  - F^"entry" ind_(b=0 " e " b'>0)
  - rho^"move" ind_(ell' != ell)
  - gamma^"conc" s_(ell')(m_y)(K_(b')/1000)^2.
$

É importante não confundir os símbolos. $R_(y,ell,theta)$, com $R$ maiúsculo, é
receita anual por MW. Já $r_y(x,a;m_y,R_y)$, com $r$ minúsculo e subscrito $y$, é a
recompensa monetária total da ação depois de combinar receita, custo e atritos. Mais
adiante, $r_R$ terá outro papel: será um resíduo adimensional que mede se a receita
foi reproduzida no ponto fixo. Portanto, $R_y$, $r_y$ e $r_R$ não são três grafias do
mesmo objeto.

Cada termo tem um sentido econômico direto. O primeiro é o lucro anual — receita por
MW $R_(y,ell',theta)$ mais o prêmio não spot $p$, menos o custo anualizado
$q_(ell',theta)$, tudo vezes a capacidade $K_(b')$. O termo $frac(chi,2)(K_(b')-
K_b)^2$ é um custo de *ajuste*: quanto mais brusca a mudança de capacidade, mais
caro, o que suaviza a expansão no tempo. $F^"entry"$ é um custo fixo de *entrada*,
pago só quando um agente sai de $b=0$ (fora do mercado) para produzir. $rho^"move"$
é um custo de *realocação*, pago quando $ell' != ell$. O último termo é um custo de
*concentração*: $s_(ell')(m_y)=sum_(theta,b)m_y(ell',theta,b)$ é a fração da
população já presente no destino, então quanto mais lotado o nó, mais caro entrar
nele.

Os parâmetros dinâmicos são $H=6$ estágios, fator de desconto $beta=0,95$,
$chi=2.000$ R\$/MW², $F^"entry"=20$ milhões de R\$, $rho^"move"=500$ milhões de R\$ e
$gamma^"conc"=5$ milhões de R\$; o prêmio $p$ herda a calibração estática e é zero no
baseline. Convém ser franco sobre o termo de concentração: ele é deliberadamente
pequeno e *não* representa uma restrição de planejamento territorial observada. Seu
papel é puramente de regularização econômica — sem ele, uma diferença mínima de
receita entre dois nós faria toda a massa despencar sobre o vencedor, um artefato de
modelo agregado que não corresponde a como projetos reais se distribuem.

== Bellman regularizado e política

Um agente que decide hoje precisa pesar o retorno imediato contra o valor de onde a
ação o deixa amanhã. O objeto que faz essa conta é o *valor da ação* $Q_y(x,a)$: a
recompensa imediata mais o valor futuro descontado do estado para onde a transição
leva,

$ Q_y(x,a)=r_y(x,a;m_y,R_y)+beta V_(y+1)(T(x,a)), $

onde $beta$ desconta o futuro e $V_(y+1)$ é o *valor do estado* no estágio seguinte.
Como $r_y$ está em reais por agente, $Q_y$ e $V_y$ também estão em reais por agente;
por isso a temperatura $tau$, que compara diferenças entre valores de ações, também
é expressa em reais.

A recursão fecha com a condição terminal

$ V_H(x)=V_6(x)=0 quad "para todo " x. $

Isso significa que não há valor de sucata, receita de continuação nem prêmio por
capacidade remanescente depois da sexta decisão. Logo, a política do estágio cinco
considera apenas sua recompensa imediata; os estágios anteriores incorporam as
recompensas futuras por meio de $beta V_(y+1)$. A condição torna o horizonte finito
bem definido, mas pode induzir efeito de borda perto do terminal. Um valor de
continuação ou horizonte maior é uma sensibilidade econômica ainda pendente.

Por que não escolher simplesmente a ação com maior $Q$? Neste problema, as ações são
discretas e várias delas podem ter valores muito próximos. Com uma melhor resposta
exata, uma variação mínima de receita ou de massa pode trocar abruptamente a ação
ótima: toda a massa de um estado iria, por exemplo, de “manter” para “expandir”, ou
de uma região para outra, mesmo que a diferença de valor fosse quase nula. Como essa
nova distribuição altera o clearing e volta a alterar as receitas, a descontinuidade
pode produzir oscilações no ponto fixo e tornar o resultado excessivamente sensível
a empates e tolerâncias numéricas.

A regularização entrópica substitui essa escolha dura por uma escolha probabilística.
Para cada estado, o agente escolhe uma distribuição $pi(.|x)$ sobre as ações e recebe,
além do valor esperado, um termo de entropia ponderado por $tau$:

$
  V_y(x)=max_(pi in Delta(cal(A)(x)))
  [sum_(a in cal(A)(x)) pi(a|x) Q_y(x,a)+tau cal(H)(pi(.|x))],
  quad cal(H)(pi)=-sum_a pi(a) log pi(a).
$

Esse problema tem uma política de valor único e contínua em $Q$: ações quase
empatadas dividem a massa gradualmente, enquanto ações claramente inferiores recebem
pouco peso. O benefício para o MFG é tornar mais suave o encadeamento “receita $arrow.r$
política $arrow.r$ distribuição $arrow.r$ clearing”, melhorando o condicionamento e a
estabilidade da iteração forward--backward. Isso não prova que o equilíbrio seja
único nem garante convergência, mas evita que pequenas perturbações provoquem saltos
puramente discretos da população.

Essa foi uma decisão metodológica durante o trabalho: a melhor resposta exata gera
um operador descontínuo justamente onde ações discretas quase empatam, enquanto o
clearing devolve receitas que mudam continuamente com a capacidade. A média suave
do `softmax` conecta esses dois objetos sem transformar diferenças numéricas mínimas
em migrações de massa de cem por cento. O benefício procurado é, portanto, duplo:
uma política única e estável para computar o ponto fixo e uma representação reduzida
de escolhas idiossincráticas entre alternativas quase equivalentes. O preço dessa
decisão é que a dispersão passa a depender de $tau$; ela não pode ser apresentada
como frequência empírica de escolhas enquanto $tau$ não for estimado ou submetido a
continuação.

A solução desse problema regularizado é o `softmax`, e seu valor é o `log-sum-exp`,
uma aproximação suave do máximo sobre as ações viáveis $cal(A)(x)$,

$ V_y(x)=tau log sum_(a in cal(A)(x)) exp(Q_y(x,a)/tau), $

e a política — a probabilidade de escolher cada ação — é o `softmax` correspondente,

$ pi_y(a|x)=frac(exp(Q_y(x,a)/tau),
  sum_(a' in cal(A)(x)) exp(Q_y(x,a')/tau)). $

O papel de $tau$ é o de uma *temperatura*: quanto maior, mais a política se espalha
entre ações de valor parecido; quanto menor, mais ela se concentra na melhor. O
baseline usa $tau=100$ milhões de R\$, um valor alto, então o objeto que calculamos
é explicitamente um *MFG entrópico regularizado*. Isso merece destaque porque muda a
interpretação: a aleatoriedade da política não é ruído acrescentado depois de otimizar
— ela *é* a solução ótima do problema regularizado. Ela também pode ser lida, de forma
reduzida, como heterogeneidade idiossincrática não observada entre projetos. Contudo,
$tau$ não foi estimado com dados de escolhas, portanto essa leitura não constitui uma
identificação empírica. Uma consequência que se deve atenção
é que um $tau$ grande suaviza escolhas de propósito, e a dispersão que veremos nos
mapas de política é, em parte, essa regularização, não necessariamente indiferença
econômica real.

Quando $tau=0$, o código substitui o `log-sum-exp` pelo máximo de Bellman e divide
probabilidade apenas entre ações numericamente empatadas:

$ pi_y(a|x)>0 arrow.r Q_y(x,a)=max_(a')Q_y(x,a'). $

Não se afirma que o equilíbrio regularizado com $tau>0$ seja idêntico ao jogo
discreto não regularizado.

#nbref([
  Notebook 07, §10, célula `[23]`, função `solve_dynamic_bellman`. Procure o ramo
  `if tau > 0`, que usa `_logsumexp` e `np.exp`, e o ramo `else`, que
  encontra `q_max`, cria `winners` e divide a probabilidade entre os índices empatados.
])

== Equação forward

A política diz o que cada agente faz; a *equação forward* agrega essas escolhas
individuais para dizer como a população inteira se move de um estágio ao seguinte. O
novo campo $hat(m)_(y+1)$ é obtido empurrando a massa atual $m_y$ através da política:
para cada estado de destino $x'$, somamos a massa que sai de cada estado $x$ e escolhe
uma ação que leva a $x'$,

$
  hat(m)_(y+1)(x')
  =sum_(x in cal(X))sum_(a in cal(A)(x))
    m_y(x) pi_y(a|x) ind_(T(x,a)=x').
$

O acento em $hat(m)$ sinaliza que é a distribuição *implicada* pela política, ainda
a comparar com a distribuição candidata $m$ que gerou as receitas e a própria
política. O chapéu não indica uma estimativa estatística: indica a saída de uma
aplicação completa da melhor resposta e da equação forward. Como a transição é
determinística e cada política soma um, a equação conserva massa por construção —
nenhum agente some ou é criado. Em memória, o notebook verifica numericamente

$ max_y abs(sum_x m_y(x)-1)=2,22 times 10^(-16), $

e, para cada estado e estágio,

$ max_(y,x) abs(sum_a pi_y(a|x)-1)=4,44 times 10^(-16). $

#nbref([
  Notebook 07, §10, célula `[23]`, função `dynamic_forward`: o laço soma
  `mass * probability` em `m_forward[stage + 1, ...]`. As verificações
  de soma da massa e da política, além da serialização detalhada, estão em §10.1,
  célula `[24]`, nos arquivos `dynamic_mfg_distribution.csv` e
  `dynamic_mfg_policy.csv`.
])

Depois da serialização em CSV e nova agregação, os limites conservadores observados
são $4,44 times 10^(-15)$ para massa e $1,11 times 10^(-15)$ para política.

== Clearing e condição de ponto fixo

Temos agora todas as peças e podemos fechar o círculo. A lógica é encadear os
operadores das subseções anteriores: em cada estágio $y$, a distribuição $m_y$ fixa a
capacidade $K_(ell,y)$; o clearing $cal(R)$ transforma capacidade em receita
$R_(y,ell,theta)$; a melhor resposta $BR$ (o Bellman regularizado) transforma receita
em política $pi_y$; e o operador forward $cal(F)$ transforma política de volta em
distribuição. O equilíbrio é o ponto fixo dessa composição,

$ m=cal(F)(BR(cal(R)(m))), $

com a igualdade valendo para a trajetória inteira $m=(m_0,dots,m_H)$ e o ponto de
partida $m_0$ fixo. Em palavras: a distribuição que sai no fim tem de ser a mesma que
entrou no começo — a população antecipa corretamente o próprio futuro.

A iteração usa

$ m^(n+1)=(1-alpha)m^n+alpha hat(m)^n, quad alpha=0,10, $

onde $n$ é o número da iteração, $m^n$ é a trajetória candidata naquela iteração,
$hat(m)^n$ é a trajetória que a melhor resposta implica a partir dela e $alpha$ é o
peso dado à nova trajetória. Com $alpha=0,10$, apenas 10% da distribuição implicada
entra na atualização e 90% da candidata anterior é preservada. Esse amortecimento
reduz oscilações, mas também pode produzir passos pequenos sem que o ponto fixo tenha
sido alcançado. Por isso o tamanho do passo $||m^(n+1)-m^n||$ não é usado sozinho
como certificado.

A certificação é feita no operador sem relaxação. O primeiro critério é

$ r_m=||hat(m)-m||_infinity
  =max_(y,ell,theta,b) abs(hat(m)_y(ell,theta,b)-m_y(ell,theta,b)), $

em que $r_m$ é o *resíduo da massa ou da distribuição*. A norma do infinito procura
a maior diferença absoluta entre $hat(m)_y(ell,theta,b)$ e $m_y(ell,theta,b)$
considerando todos os estágios, localizações, tipos e bins. Em palavras, ele localiza
a célula em que a população implicada pela política mais diverge da população que os
agentes anteciparam. Como $m$ é uma fração, $r_m$ também é uma fração da população.
Se $r_m=0$, as duas distribuições coincidem em todas as células; se $r_m$ é pequeno,
nenhuma célula diverge muito.

O segundo critério é

$ r_C=max_(y,ell) frac(abs(hat(K)_(y,ell)-K_(y,ell)),K_(y,ell)+1), $

em que $r_C$ é o *resíduo da capacidade*. Apesar do subscrito $C$, a capacidade é
denotada por $K$. Primeiro convertemos $m$ em $K_(y,ell)$ MW e $hat(m)$ em
$hat(K)_(y,ell)$ MW; depois calculamos a diferença relativa entre essas duas
capacidades e tomamos o pior caso entre todos os estágios e localizações. O $+1$ no
denominador impede divisão por zero quando a capacidade é nula ou muito pequena.
Assim, $r_C$ responde a uma pergunta diferente de $r_m$: mesmo uma diferença
populacional pequena pode ser material em MW se estiver concentrada em bins grandes.

O terceiro critério é

$ r_R=max_(y,ell,theta)
  frac(abs(hat(R)_(y,ell,theta)-R_(y,ell,theta)),
       abs(R_(y,ell,theta))+1). $

Aqui $r_R$ é o *resíduo da receita*. A receita $R_(y,ell,theta)$ é calculada pelo
clearing usando a capacidade proveniente de $m$. Em seguida, refazemos o clearing
com a capacidade proveniente de $hat(m)$ e obtemos $hat(R)_(y,ell,theta)$. O resíduo
compara as duas receitas em termos relativos e retém a maior discrepância entre
estágios, localizações e tipos solares; novamente, o $+1$ evita divisão por zero.
Portanto, $r_R$ verifica se a realimentação econômica também fechou: não basta a
distribuição parecer próxima se essa pequena diferença alterar muito preços,
curtailment e receita por MW.

Os três resíduos seguem a mesma regra: quanto menores, mais perto do ponto fixo;
zero seria consistência exata. Eles não são probabilidades de erro nem intervalos de
confiança. São distâncias determinísticas entre a trajetória candidata e a trajetória
produzida pelo próprio operador. $r_C$ e $r_R$ são adimensionais por construção;
$r_m$ é uma fração da população e também não possui unidade física.

As tolerâncias são $2 times 10^(-4)$ para $r_m$ e
$5 times 10^(-3)$ para $r_C$ e $r_R$. A trajetória só é certificada quando os três
critérios passam simultaneamente: consistência da população, da capacidade física e
da receita econômica.

#nbref([
  Notebook 07, §10, célula `[23]`, função `solve_dynamic_mfg`. Procure os blocos
  finais `residual_m`, `residual_C`, `residual_R` e
  `certified = bool(residual_m < DYN_TOL_M and residual_C < DYN_TOL_C and
  residual_R < DYN_TOL_R)`. A série por iteração é gravada em §10.1, célula
  `[24]`, como `dynamic_mfg_convergence.csv`.
])

A pergunta que a próxima figura precisa responder é dupla: o ponto fixo dinâmico de
fato existe (a iteração converge nas tolerâncias) e, existindo, para onde ele leva a
capacidade ao longo dos estágios. É o resultado central da seção, e vale lê-lo com o
mesmo ceticismo aplicado ao benchmark estático — convergência é consistência interna,
não veracidade externa.

A pergunta que a próxima figura precisa responder é dupla: o ponto fixo dinâmico de
fato existe (a iteração converge nas tolerâncias) e, existindo, para onde ele leva a
capacidade ao longo dos estágios. É o resultado central da seção, e vale lê-lo com o
mesmo ceticismo aplicado ao benchmark estático — convergência é consistência interna,
não veracidade externa.

#imgfig(
  [Certificação do MFG dinâmico e trajetória de capacidade.],
  "outputs/figures/rel_dynamic_mfg_fixed_point.png",
  note: [
    À esquerda, em escala logarítmica, os resíduos do operador (antes da relaxação) e
    o passo numérico caindo com as iterações; à direita, a capacidade $K_(ell,y)$ de
    cada subsistema ao longo dos estágios. No plano numérico, a iteração 67 encerra
    com $r_m=2,00 times 10^(-4)$, $r_C=2,19 times 10^(-4)$ e
    $r_R=2,61 times 10^(-6)$, os três dentro das tolerâncias. Existe, portanto, uma
    trajetória internamente consistente nos três níveis, não apenas uma sequência de
    passos pequenos.

    O painel direito mostra *capacidade instalada*, não geração. Com a ligação
    NE--SE, os preços-sombra médios iniciais ficam próximos no SE, S e NE
    (274,97; 275,34; 274,96 R\$/MWh), enquanto a produtividade solar faz a receita
    `high` do NE atingir R\$ 764,16 mil/MW/ano, acima dos R\$ 685,73 mil do SE.
    Assim, SE e NE crescem em conjunto: o SE vai de 19,077 para 26,345 GW e o NE de
    11,789 para 20,847 GW. A massa terminal fica em 48,43% no SE, 36,94% no NE,
    9,32% no S e 5,31% no N.

    Duas ressalvas críticas equilibram o resultado. A primeira: 53,685 GW fica
    *abaixo* dos 80,764 GW do planejador social e *acima* dos 38,954 GW do privado
    estático, o que faz sentido — o agente dinâmico enxerga valor futuro que o privado
    míope não vê, mas continua sem internalizar as externalidades que só o planejador
    considera. A segunda, e é a que mais importa: nada disso valida os parâmetros. O
    número terminal é fortemente condicionado por $tau$, pelos custos de entrada e
    realocação e pelo horizonte de seis estágios. A vantagem de receita do Sudeste/CO
    e do Nordeste também herda a rede reduzida, os perfis e a ausência de contratos;
    mudar esses objetos pode alterar a hierarquia locacional. O que a figura certifica
    é o ponto fixo desta formulação. Os dados estão em `dynamic_mfg_path.csv` e
    `dynamic_mfg_convergence.csv`, gerados no notebook 07, §10--10.1, células
    `[23]`--`[24]`; a figura é reproduzida no notebook 08, §12, célula `[13]`.
  ]
)

== Trajetória agregada

#text(size: 8.35pt)[
#table(
  columns: (0.72fr, 1.02fr, 1.02fr, 1.02fr, 1.02fr),
  inset: 3.8pt,
  align: right,
  [Subs.], [Massa inicial], [Massa terminal], [Cap. inicial GW], [Cap. terminal GW],
  [SE], [48,83%], [48,43%], [19,077], [26,345],
  [S], [15,15%], [9,32%], [5,919], [4,030],
  [NE], [30,18%], [36,94%], [11,789], [20,847],
  [N], [5,85%], [5,31%], [2,284], [2,464],
  [*Total*], [*100%*], [*100%*], [*39,069*], [*53,685*],
)
]

A capacidade pode aumentar mesmo quando a massa de uma região cai, pois os agentes
remanescentes podem migrar para bins maiores. Esse é o motivo para reportar massa e
capacidade separadamente.

A trajetória agregada esconde *como* a expansão acontece — se todos crescem um pouco
ou se alguns saltam para bins altos. A figura seguinte abre a distribuição inteira
para mostrar esse mecanismo, e o motivo de olhá-la é justamente distinguir crescimento
por adesão de crescimento por intensificação.

#imgfig(
  [Trajetória completa da distribuição por localização e capacidade individual.],
  "outputs/figures/rel_dynamic_mfg_distribution_trajectory.png",
  note: [
    Esta figura também não mostra geração horária: mostra onde a *massa de agentes*
    está na grade de capacidade. Cada painel é um subsistema; o eixo horizontal é o
    estágio $y=0,dots,6$, o vertical é a capacidade individual $K_b$ e a cor é a massa $sum_theta m_y(ell,
    theta,b)$ naquele nó e bin, com escala comum aos quatro. Toda a população começa
    concentrada em 372 MW por agente e, ao longo do horizonte, SE e NE empurram massa
    para bins maiores — é intensificação, não apenas entrada de mais agentes.

    No bin inicial de 372 MW, o ajuste esperado no estágio zero é +57,85 MW (`high`)
    e +25,65 MW (`low`) no SE, contra +64,84 e +36,32 MW no NE. No Sul, os valores
    são +46,21 e +5,39 MW; no Norte, +47,72 e +7,81 MW. O Nordeste apresenta,
    portanto, o sinal de expansão mais forte no início, coerente com sua
    receita solar por MW. O ponto que exige cuidado de leitura
    é uma armadilha visual real: uma célula escura pode significar duas coisas opostas
    — pouca massa naquela região, ou massa espalhada por muitos bins, cada um com
    pouca cor. As duas se parecem no mapa e não devem ser confundidas; por isso esta
    figura só se interpreta corretamente ao lado da massa locacional a seguir, que
    resolve a ambiguidade somando sobre os bins.
  ]
)

Esta figura é o contraponto que resolve a ambiguidade da anterior: em vez de bins,
mostra a fração de agentes por região, respondendo diretamente à pergunta "para onde
a população migra".

#imgfig(
  [Evolução da fração da população por subsistema.],
  "outputs/figures/rel_dynamic_mfg_location_mass.png",
  note: [
    A altura empilhada é sempre um e cada faixa é a fração da população em uma
    localização, integrando tipos e bins. O movimento é claro: o Nordeste ganha
    massa relativa, de 30,18% para 36,94%; o SE fica quase estável (48,83% para
    48,43%); Sul e Norte recuam. No bin inicial, a probabilidade de realocação é
    baixa no SE (2,26% `high`; 1,90% `low`) e no NE (0,26%; 0,34%), mas elevada
    no Sul (27,53%; 16,22%) e no Norte (16,32%; 12,20%). Assim, SE e NE retêm
    agentes, enquanto as regiões com margens menores perdem participação.

    O fato de o Nordeste oferecer margem maior não faz todos os agentes escolherem o
    mesmo caminho. Primeiro, o tipo `high` produz mais e tem custo diferente do tipo
    `low`, de modo que dois projetos no mesmo estado físico podem avaliar a expansão
    de forma distinta. Segundo, mudar de região custa R\$ 500 milhões e ajustar
    capacidade incorre em uma penalidade quadrática calibrada em $2.000$ R\$/MW²;
    localização e porte atuais importam, e a decisão ótima depende do ponto de partida.
    Terceiro, quando muitos agentes ocupam
    um destino, a penalidade de concentração cresce com a massa local e com o quadrado
    da capacidade. Além disso, a própria expansão reduz a receita solar local por
    canibalização. O destino que parece melhor para um agente isolado deixa de ser
    universalmente melhor quando a população se move.

    A regularização entrópica com $tau>0$ também transforma a escolha em uma política
    probabilística: ações com valores próximos recebem massa positiva, em vez de toda
    a população executar um único `argmax`. Esse mecanismo é estilizado, mas traduz
    aspectos reais — projetos têm produtividades, custos, terrenos, pontos de conexão,
    contratos e expectativas diferentes, e decisões de grande porte não se deslocam
    sem custo. Assim, a dispersão não é um acidente do gráfico: ela nasce da
    heterogeneidade permanente, dos custos de ajuste e mudança, da congestão de
    destinos e da realimentação de mercado.

    O custo de mudança de R\$ 500 milhões impede que essa migração seja automática: a
    vantagem do destino precisa compensar o custo presente por meio da soma descontada
    das margens futuras. A política entrópica ainda distribui alguma probabilidade por
    ações não modais, de modo que esses números não significam que um quarto das usinas
    físicas do Nordeste seja transportado de região; significam que a massa de
    decisões de investimento se redistribui no modelo. Cruzando com a tabela agregada,
    aparece a distinção que justifica reportar massa e capacidade separadamente: o
    Nordeste *ganha* participação e capacidade (11,789 para 20,847 GW), enquanto a
    capacidade do SE também cresce mesmo com massa quase estável. Um alerta de
    interpretação para não ler
    demais na figura: isto é realocação de *agentes de investimento*, não fluxo físico
    de energia entre subsistemas — a rede da Seção 4 é que move elétrons, esta figura
    move decisões.
  ]
)

== Política de localização e capacidade

As figuras anteriores mostraram o *resultado* do movimento; as duas próximas mostram
a *regra de decisão* que o produz — a política em si. Queremos ver de onde vem a
migração para o corredor SE--NE: retenção nesses nós ou saída das outras regiões? E
onde a política decide crescer versus encolher?

#imgfig(
  [Política no estágio inicial: realocação e ajuste esperado.],
  "outputs/figures/rel_dynamic_mfg_policy.png",
  note: [
    As colunas separam os tipos `high` e `low`; na linha de cima a cor é a
    probabilidade de mudar de localização ($ell' != ell$) e na de baixo é o ajuste
    esperado de capacidade $EE[Delta K|x]$ em MW, vermelho para expansão e azul para
    redução. Esse ajuste esperado é a média ponderada das mudanças possíveis,
    $EE[Delta K|x]=sum_a pi_y(a|x)(93 delta)$ MW, pois cada ação altera o bin em
    $delta=-1,0,+1$ e cada bin vale 93 MW. As localizações atuais estão nas linhas e
    os bins atuais nas colunas. A leitura econômica é coerente com a migração
    observada: estados do SE e do NE quase não querem sair e ainda expandem em quase
    toda a grade, enquanto Sul e Norte mostram maior propensão a migrar.

    A própria fórmula da recompensa explica por que a realocação cresce com o bin. A
    vantagem corrente de trocar uma origem $ell$ por um destino $ell'$ contém a diferença
    de margens multiplicada por $K_b$,
    $[(R_(ell')-q_(ell'))-(R_ell-q_ell)]K_b$, enquanto o custo de mudança
    $rho^"move"$ é fixo. Quanto maior o produtor, mais a diferença de margem pesa e
    menor fica, quando dividido pelos MW do projeto, o peso do custo fixo de R\$ 500
    milhões. No sentido inverso,
    sair de SE ou NE significa abandonar margens altas e ainda pagar esse custo,
    daí as probabilidades pequenas nessas duas linhas.

    As cores também mostram a realimentação de canibalização ao longo do horizonte.
    A expansão acumulada derruba a receita local e reduz gradualmente o incentivo a
    continuar crescendo, razão pela qual as trajetórias perdem inclinação. O movimento
    combina retenção no corredor SE--NE, saída seletiva do Sul e do Norte e freio por
    canibalização. Há um artefato de borda que vale
    apontar para não ser mal lido: no bin máximo o ajuste esperado fica negativo
    simplesmente porque a ação $+1$ é inviável e a política ainda coloca alguma
    probabilidade em reduzir — é a fronteira da grade, não um sinal de desinvestimento
    econômico. E a cautela de sempre com este modelo: a regularização entrópica com
    $tau$ alto borra as escolhas, então as cores intermediárias são em boa parte a
    temperatura da política, não heterogeneidade observada nem ruído numérico.
  ]
)

Por fim, uma visão individual da política. É a mais intuitiva das figuras, mas também
a mais fácil de sobreinterpretar, e vale usá-la sabendo exatamente o que ela pode e
o que não pode dizer.

#imgfig(
  [Trajetórias de ações modais para agentes representativos.],
  "outputs/figures/rel_dynamic_mfg_modal_paths.png",
  note: [
    Cada linha fixa uma localização inicial, cada coluna é um estágio, e cada célula
    mostra o destino e a capacidade da *ação modal* — a mais provável — com os painéis
    separando tipos. Como narrativa, ela é útil: o `high` modal sobe um bin em todos
    os estágios e regiões, chegando a 930 MW. O `low` termina em 837 MW no SE e no
    NE e permanece em 372 MW no Sul e no Norte. Esse contraste é o
    resultado esperado da margem: o tipo `high` produz mais e custa menos, e o
    corredor SE--NE sustenta receita por MW maior, de modo que suporta expansão por mais
    estágios antes de o custo de ajuste e a canibalização dominarem.

    Há uma aparente contradição instrutiva: nenhuma trajetória modal troca de região,
    embora a figura de massa mostre migração agregada para o Nordeste. Não é erro. Em
    cada estado, “ficar” pode ser a ação individual mais provável e ainda assim a
    *soma* das probabilidades de várias ações de realocação mover massa relevante. No
    Sul `high` inicial, por exemplo, ficar e expandir é a ação modal com 41,04%,
    mas todas as ações de mudança somam 27,53%. A equação forward transporta essa
    massa; a trajetória modal a descarta.

    Por isso o limite é fundamental e precisa ser dito sem meias palavras: a ação
    modal *joga fora* toda a probabilidade das ações
    não modais, e somar essas trajetórias representativas *não* reconstrói a
    distribuição $m_y$ — uma população em que 51% sobe e 49% desce apareceria aqui
    como "todos sobem", o que é falso no agregado. Qualquer conclusão sobre a
    população tem de vir da equação forward e dos mapas de distribuição; esta figura é
    ilustração pedagógica, não evidência quantitativa.
  ]
)

#resultbox([Resultado central do MFG dinâmico])[ 
  Existe uma trajetória internamente consistente do MFG entrópico regularizado nas
  tolerâncias declaradas. Ela desloca massa para o Nordeste, aumenta o bin médio
  de capacidade e leva o total de 39,069 a 53,685 GW. Localização, capacidade e valor
  futuro são escolhidos em conjunto. O rótulo “certificado” refere-se ao ponto fixo do modelo,
  não à validade externa dessa trajetória.
]

= Síntese dos resultados

#text(size: 8.65pt)[
#table(
  columns: (1.55fr, 1.25fr, 2.55fr),
  inset: 4pt,
  align: left,
  [Objeto], [Resultado], [Leitura correta],
  [Painel], [35.040 linhas], [8.760 horas completas por subsistema em 2025.],
  [Captura solar], [0,640 a 0,729], [CMO ponderado pela solar abaixo da média simples nos quatro subsistemas.],
  [Reservatório], [$4,39 times 10^(-11)$ MWh], [erro máximo da equação de movimento; fechamento cíclico exato.],
  [Cenário seco], [52,93 TWh de déficit], [choque uniforme de -30% na afluência, demais entradas fixas.],
  [Benchmark estático], [38,954 GW], [aproxima o alvo de 39,069 GW e é certificado por três resíduos.],
  [Total fixo], [39,069 GW], [contrafactual de realocação; não é equilíbrio privado.],
  [Social limitado], [80,764 GW], [nenhum teto regional ativo; resultado ainda condicional.],
  [MFG dinâmico], [53,685 GW terminal], [MFG entrópico regularizado certificado na iteração 67.],
  [Conservação], [$4,44 times 10^(-15)$], [erro máximo após serializar e reagrupar a distribuição; em memória, $2,22 times 10^(-16)$.],
)
]

Os resultados contam uma história coerente. A solar possui fator de captura abaixo
de um; o clearing mostra que a resposta do preço depende de água, déficit e rede;
o laboratório finito revela um limiar de adequação; o benchmark estático fecha o
nível calibrado; o planejador e o MFG alteram esse nível porque internalizam objetos
diferentes e, no caso dinâmico, valor futuro e realocação.

= Limitações, identificação e extensões

== Limitações empíricas

- A capacidade `p99` mede produção efetiva extrema, não placa instalada nem
  disponibilidade contratual.
- O CMO observado sustenta diagnósticos de custo, mas não substitui PLD, receita
  contratual ou fluxo de caixa de um projeto.
- Os quatro medoides comprimem variabilidade diária e não preservam sequências de
  seca, enchimento ou rampas entre dias.
- A rede tem cinco nós contábeis e cinco interfaces fixas; perdas, contingências e
  restrições internas não estão representadas.

== Limitações do clearing

- Hidrologia é agregada em energia equivalente; não há cascatas, curvas de queda,
  evaporação, defluência mínima ou coordenação plurianual.
- Custos térmicos são blocos reduzidos ancorados em estatísticas de CMO, não ofertas
  observadas de cada usina.
- Curtailment usa epsilon de desempate. Compensação financeira depende de causa e
  contrato e precisa de modelagem própria.
- O regularizador quadrático melhora unicidade e estabilidade, mas entra no preço
  marginal. Sensibilidades a $epsilon_x$ devem acompanhar qualquer exercício
  normativo.

== Limitações de investimento e MFG

- O custo locacional e a produtividade dos dois tipos são proxies. Embora nenhum teto
  regional esteja ativo no baseline atual, a fronteira de capacidade continua sendo
  uma hipótese que requer sensibilidade.
- O prêmio não spot é um único escalar. Contratos devem ser representados por
  volume, preço, prazo, perfil, indexação e penalidades.
- O horizonte dinâmico possui seis estágios abstratos. Ele organiza decisões
  sequenciais, mas não está ainda calibrado como seis anos civis específicos.
- $tau=100$ milhões de R\$ suaviza a política. Uma análise robusta deve traçar a
  continuação em $tau$ e comparar com o caso $tau=0$.
- O certificado verifica o ponto fixo encontrado, mas não demonstra unicidade
  global. Inicializações alternativas e análise de estabilidade continuam
  necessárias.

== Extensões prioritárias

#enum[
  substituir as proxies p99 por cadastro de capacidade, disponibilidade e conexão;
][
  incorporar PLD observado e contratos sem confundi-los com CMO ou $lambda$;
][
  usar cenários hidrológicos cronológicos com reservatórios acoplados entre dias;
][
  estimar limites horários de intercâmbio e, quando possível, uma rede nodal mais
  detalhada;
][
  executar sensibilidades do planejador aos tetos regionais, VOLL, CAPEX e valor da
  água;
][
  estudar continuação em $tau$, múltiplas inicializações e seleção de equilíbrio;
][
  validar previsões locacionais fora da amostra, separando calibração do total e
  avaliação regional.
]

= Conclusão

O relatório organiza uma cadeia de objetos verificáveis. Carga residual e captura
explicam o mecanismo; o clearing reúne as restrições físicas em um programa
quadrático convexo; as KKT dão interpretação ao preço-sombra e ao valor da água; o
laboratório finito mostra como entrada altera restrições ativas; e o MFG fecha
Bellman, distribuição forward e clearing em uma trajetória consistente.

Os elementos que sustentam essa cadeia são curtailment explícito, estoque hídrico
com fechamento cíclico, sinal do dual determinado pela convenção do balanço,
distinção entre passo relaxado e resíduo do operador, representação baricêntrica sem
falsa mistura discreta, planejador independente e nomenclatura explícita do MFG
entrópico. A rastreabilidade das fontes, a classificação dos parâmetros e a
separação entre externalidade e diferença de especificação delimitam o que os
resultados permitem concluir.

#resultbox([Mensagem final])[ 
  A conclusão defensável não é um número isolado de capacidade. É a existência de
  uma cadeia auditável: dados horários determinam perfis; perfis e capacidade
  determinam o clearing; o clearing determina receitas e preços-sombra; agentes
  respondem a esses sinais; e a distribuição resultante deve reproduzir o campo
  usado para calculá-los. Quando essa última igualdade é afirmada, ela vem
  acompanhada do resíduo que a certifica.
]

#pagebreak(weak: true)
