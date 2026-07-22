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

#align(center)[
  #text(17.2pt, weight: "bold")[MFG - Renewable Energy]

  #v(0.24em)
  #text(13.0pt, weight: "semibold", fill: rgb("2f5f7f"))[
    Relatório técnico - versão 2
  ]

  #v(0.52em)
  #text(10.1pt)[João Felipe Vilas Boas]

  #v(0.18em)
  #text(8.8pt, fill: gray.darken(20%))[
    Modelo, resultados e auditoria computacional dos notebooks 07 e 08
  ]
]

#v(0.80em)

#block(inset: 9pt, stroke: 0.65pt + gray, radius: 5pt, fill: luma(250))[
  *Resumo.* Este relatório desenvolve uma estrutura empírico-matemática para
  estudar como a expansão de fontes renováveis intermitentes altera a carga
  residual, o despacho hidrotérmico, o preço-sombra locacional, o valor da água e
  os incentivos de investimento no Sistema Interligado Nacional. A versão 2
  preserva o núcleo do relatório original, mas substitui resultados preliminares
  por uma implementação integrada e auditável. O painel contém 35.040 observações
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

#pagebreak()

#align(center)[#text(14pt, weight: "bold")[Sumário]]
#v(0.35em)
#columns(2, gutter: 1.0cm)[
  #set text(size: 8.35pt)
  #outline(title: none)
]

#pagebreak()

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

A pergunta principal permanece a do relatório original:

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

== O que é preservado e o que é corrigido

O cerne conceitual do `report.typ` está presente nos notebooks 07 e 08. A
equivalência não significa repetir números antigos: significa preservar as
relações matemáticas e econômicas, corrigindo formulações que eram apenas
preliminares.

#text(size: 8.7pt)[
#table(
  columns: (1.45fr, 1.15fr, 2.75fr),
  inset: 4pt,
  align: left,
  [Núcleo do relatório original], [Situação na versão 2], [Implementação e correção],
  [Carga residual, curva do pato e captura], [Preservado],
  [Painel anual de 2025 e figuras do notebook 08; CMO permanece separado de PLD e de $lambda$.],
  [Hidrologia proxy e valor da água], [Fortalecido],
  [Reservatório com estado $S_t$, afluência, turbinamento, vertimento e $S_24=S_0$; a passagem física para energia equivalente é explícita.],
  [Despacho, Lagrangiana e KKT], [Corrigido],
  [Um único primal é compartilhado pelas camadas; o sinal do dual segue o balanço $L-g=0$ e o curtailment pertence ao problema.],
  [Dois produtores solares], [Preservado],
  [Incumbente e entrante usam o mesmo clearing, com receita total, receita capturada, lucro, déficit e corte separados.],
  [Energia diária versus potência horária], [Preservado e explicitado],
  [A desigualdade agregada não é usada como prova de adequação em cada hora.],
  [Seco, base e úmido], [Recalculado],
  [Apenas a afluência muda; demanda, capacidade, custo e rede ficam fixos.],
  [Orçamento fixo, entrada e ótimo social], [Reformulado],
  [O total fixo é um contrafactual de realocação, não um falso equilíbrio. O privado é calibrado e o planejador é resolvido independentemente.],
  [MFG locacional], [Substancialmente fortalecido],
  [Estado, ação, política, Bellman, forward, clearing e distribuição formam um ponto fixo dinâmico certificado por resíduos não relaxados.],
)
]

#resultbox([Conclusão da auditoria conceitual])[ 
  Após as explicitações acima, não há seção matemática central do relatório antigo
  sem contrapartida nos notebooks 07 e 08. O notebook 06 foi usado como referência
  histórica da implementação, mas o 07 substitui seu orçamento hídrico diário por
  uma dinâmica de estoque, corrige o certificado numérico e acrescenta o MFG
  dinâmico completo. O notebook 08 transforma esses objetos em figuras, sem
  reestimar ou alterar o modelo.
]

== Como ler resultados e certificados

Há três níveis de afirmação que não devem ser confundidos:

- um *fato descritivo* é calculado diretamente do painel, como correlação ou fator
  de captura;
- um *resultado de modelo* depende de parâmetros e proxies, como a capacidade do
  planejador;
- um *certificado numérico* apenas afirma que a solução satisfaz o operador e as
  tolerâncias declaradas; ele não transforma a calibração em evidência causal.

Da mesma forma, CMO, PLD e preço-sombra são objetos diferentes. Usamos
$CMO^"obs"$ para os diagnósticos empíricos; a base atual não contém PLD observado;
e $lambda^"model"$ é o dual do balanço no clearing acadêmico.

= Dados, notação e compressão temporal

== Índices, estados e unidades

O tempo operacional é horário, $t=0,dots,23$. Os dias representativos são
$d in cal(D)$, os subsistemas são

$ cal(L) = ("SE", "S", "NE", "N"), $

os blocos térmicos são $b in cal(B)^"th"$ e os tipos solares são
$theta in ("high", "low")$. No MFG dinâmico, $y=0,dots,H$ denota o estágio de
decisão e não a hora do despacho.

#keybox([Unidades e objetos que exigem cuidado])[ 
  Capacidade e potência são medidas em MW; energia horária em MWh; estoques e
  afluências em MWh-equivalentes; receitas anuais por unidade de capacidade em
  R\$/MW/ano; e $lambda$ e CMO em R\$/MWh. A capacidade `p99` é uma proxy
  operacional baseada no percentil 99 da geração observada, não capacidade nominal
  cadastrada. Os pesos dos dias representativos são reescalados para 365 dias e,
  portanto, anualizam o modelo.
]

== Painel horário de 2025

O notebook 07 lê o painel construído pelo notebook 00, padroniza as colunas e faz
uma junção à esquerda do CMO por `(instante, subsistema)`. A única lacuna diária do
cache de preço é preenchida dentro de cada subsistema por propagação temporal; não
há propagação entre regiões. O painel consolidado possui 8.760 horas por
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
cadastro físico, são:

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

Resolver 365 despachos de 24 horas a cada iteração seria desnecessariamente caro.
O número de clusters é escolhido entre 4 e 12 combinando silhueta, estação e
regime de calendário. Cada cluster é representado pelo dia real mais próximo do
centro, isto é, por um medoide; nenhum perfil horário sintético é inventado.

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

#imgfig(
  [Dias representativos e pesos anualizados.],
  "outputs/figures/rel_dias_representativos.png",
  note: [
    *Como ler.* Cada barra é um medoide real; sua altura é o número de dias do ano
    representado pelo cluster. Cor identifica estação e a anotação identifica o
    regime de calendário. *Resultado.* Os pesos 71, 84, 104 e 106 somam 365, de
    modo que receitas e custos operacionais podem ser anualizados sem multiplicar
    novamente por 365. *Limite.* Medoides preservam perfis diários observados, mas
    não preservam uma cronologia sazonal entre dias; por isso cada reservatório é
    fechado ciclicamente dentro do bloco.
  ]
)

#pagebreak(weak: true)

= Mecanismo empírico: carga residual, preço e captura

== Carga líquida, carga residual e rampa

Se $L_(ell,t)$ é a carga do painel, $G_(ell,t)^s$ a solar observada e
$G_(ell,t)^w$ a eólica observada, definimos

$ L_(ell,t)^"net" = L_(ell,t) - G_(ell,t)^s, $

$ L_(ell,t)^"res" = L_(ell,t) - G_(ell,t)^s - G_(ell,t)^w, $

e a rampa residual

$ Delta L_(ell,t)^"res" = L_(ell,t+1)^"res" - L_(ell,t)^"res". $

A carga residual é o objeto que aproxima quanto o restante do sistema precisa
fornecer depois das fontes variáveis. Ela não é uma demanda econômica e não inclui,
por si só, restrições de rede, rampa térmica ou estoque hídrico.

#imgfig(
  [CMO observado versus carga residual por subsistema.],
  "outputs/figures/rel_cmo_vs_carga_residual.png",
  note: [
    *Como ler.* Cada ponto é uma hora; o eixo horizontal mostra carga menos solar
    menos eólica, e o vertical mostra o CMO observado. A reta é uma regressão
    descritiva e o número em cada painel é a correlação de Pearson. *Resultado.* As
    correlações são 0,581 no Norte, 0,412 no Sudeste/CO e 0,281 no Sul. No Nordeste
    ela é -0,122: a oferta eólica local frequentemente excede a carga regional, de
    modo que intercâmbio, excedente renovável e restrições do restante do sistema
    rompem a relação local simples. *Cuidado.* Correlação não identifica a curva de
    oferta nem prova causalidade; ela apenas documenta o canal empírico que o
    clearing tornará estrutural.
  ]
)

#imgfig(
  [Perfil horário do CMO observado com média e dispersão.],
  "outputs/figures/rel_cmo_perfil_horario_erro.png",
  note: [
    *Como ler.* O ponto central é a média do CMO em cada hora do dia; a barra é mais
    ou menos um desvio-padrão, calculado separadamente por subsistema. *Resultado.*
    A dispersão é grande em torno da média, especialmente nas horas críticas. Isso
    mostra por que uma única curva média não representa todos os regimes
    operacionais. *Cuidado.* As barras não são intervalos de confiança e não
    implicam normalidade; medem variabilidade histórica dentro de 2025.
  ]
)

#imgfig(
  [Curva do pato no Sudeste/CO: carga, líquida de solar e residual.],
  "outputs/figures/rel_curva_pato.png",
  note: [
    *Como ler.* A curva preta é a carga média; a laranja subtrai solar; a verde
    também subtrai eólica. O vale residual médio ocorre por volta das 12h, em
    aproximadamente 30,1 GW. *Resultado.* Quando o sol se põe, a maior subida média
    de uma hora é cerca de 5,17 GW entre as horas 16 e 17. A solar reduz energia
    diurna, mas aumenta a necessidade de flexibilidade no entardecer. *Cuidado.* É
    um perfil médio; máximos em dias individuais podem ser maiores e são tratados
    pelos medoides no clearing.
  ]
)

== Receita capturada e canibalização

Para uma fonte com geração $G_t^r$ e uma série de preço $P_t$, o preço capturado é

$ P^"cap,r" = frac(sum_t P_t G_t^r, sum_t G_t^r), $

e o fator de captura é

$ F^"cap,r" = frac(P^"cap,r", T^(-1) sum_t P_t). $

Valores abaixo de um significam que a fonte gera mais em horas de preço abaixo da
média. No diagnóstico empírico usamos $P_t=CMO_t^"obs"$; isso não transforma CMO
em PLD.

O mecanismo estrutural pode ser visto com

$ lambda_t = Lambda(L_t-K^"agg" a_t^s-G_t^w), quad Lambda' > 0. $

Para um projeto pequeno com perfil $a_(i,t)^s$,

$ R_i^"cap" = frac(sum_t lambda_t a_(i,t)^s, sum_t a_(i,t)^s), $

e, mantendo os perfis fixos,

$ frac(partial R_i^"cap", partial K^"agg")
  = - frac(sum_t a_(i,t)^s a_t^s Lambda'_t, sum_t a_(i,t)^s) <= 0. $

A desigualdade é estrita quando há coincidência entre o projeto e a solar agregada
em horas nas quais o preço responde à carga residual.

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

#imgfig(
  [Fatores de captura solar e eólico usando CMO observado.],
  "outputs/figures/rel_fatores_captura.png",
  note: [
    *Como ler.* A linha horizontal em um separa preço capturado igual à média de
    preço capturado abaixo ou acima da média. *Resultado.* A solar fica abaixo de
    um nos quatro subsistemas, entre 0,640 e 0,729; a eólica fica acima de um, entre
    1,044 e 1,309. O resultado é compatível com coincidência solar mais forte e
    complementaridade temporal da eólica na amostra. *Cuidado.* Trata-se de fator
    calculado com CMO, não de retorno financeiro realizado; contratos, PLD,
    encargos e diferenças entre geração centralizada e distribuída não estão nessa
    conta.
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

Todas as camadas de decisão chamam o mesmo operador físico. Dada a capacidade
solar por localização e tipo, o clearing realiza a transformação

$ K arrow.r overline(g)^s(K) arrow.r
  (g^s,c^s,g^w,c^w,h,S,v,n,F,u,lambda). $

Assim, diferenças entre o laboratório finito, o benchmark privado, o
contrafactual fixo, o planejador e o MFG não vêm de trocar silenciosamente o
despacho. Elas vêm apenas de quem escolhe capacidade e de como a decisão retorna ao
campo agregado.

== Rede e variáveis de decisão

Há quatro nós físicos e um hub contábil $I$. As arestas são

$ I arrow.l.r "SE", quad I arrow.l.r "NE", quad I arrow.l.r "N",
  quad "SE" arrow.l.r "S". $

O hub não possui carga nem geração e satisfaz conservação algébrica dos fluxos. Os
limites reduzidos são 8,0 GW em I-SE, 8,0 GW em I-NE, 4,5 GW em I-N
e 6,0 GW em SE-S. Não existem ligações diretas N-NE ou N-SE no
grafo reduzido.

#imgfig(
  [Topologia reduzida de intercâmbio utilizada pelo clearing.],
  "outputs/figures/rel_topologia_intercambio.png",
  note: [
    *Como ler.* Os quatro círculos coloridos são subsistemas com carga e geração;
    o círculo cinza $I$ é apenas um nó de contabilidade. Cada seta dupla representa
    uma única variável de fluxo com sinal e limite bilateral, não dois fluxos
    independentes. *Resultado.* Norte, Nordeste e Sudeste/CO trocam energia via
    hub; o Sul se conecta ao Sudeste/CO. *Cuidado.* Esta rede é deliberadamente
    agregada. Ela representa congestionamento inter-regional, mas não perdas,
    restrições internas, contingências ou limites horários observados do ONS.
  ]
)

Para cada $(d,ell,t)$, as variáveis do primal são:

- $g^s$ e $c^s$: solar entregue e solar cortada;
- $g^w$ e $c^w$: eólica entregue e eólica cortada;
- $h$, $S$ e $v$: geração hídrica, estoque e vertimento;
- $n_b$: geração do bloco térmico $b$;
- $u$: déficit, isto é, energia não atendida penalizada;
- $F_e$: fluxo orientado na aresta $e$.

== Oferta renovável e curtailment explícito

Para a solar, a disponibilidade depende da capacidade endógena:

$ overline(g)_(d,ell,t)^s(K)
  = a_(d,ell,t)^s sum_theta xi_theta K_(ell,theta), $

onde $xi_"high"=1,08$ e $xi_"low"=0,94$. O despacho separa energia entregue e
cortada:

$ g_(d,ell,t)^s + c_(d,ell,t)^s = overline(g)_(d,ell,t)^s(K), $

$ g_(d,ell,t)^w + c_(d,ell,t)^w = overline(g)_(d,ell,t)^w. $

Essa igualdade é conceitualmente importante. Curtailment não é calculado depois
como uma diferença potencialmente negativa; é uma variável não negativa do
problema e participa das condições de complementaridade.

== Da hidráulica física ao reservatório em energia equivalente

Em uma usina $j$, a potência física é

$ P_(j,t)^h = eta_(j,t) rho^w g H_(j,t)^"net" Q_(j,t)^"tur", $

com altura líquida

$ H_(j,t)^"net" = Z_j^"mon"(V_(j,t))
  - Z_j^"jus"(Q_(j,t)^"tur"+Q_(j,t)^"spill")
  - H_(j,t)^"loss". $

Uma representação por usina exigiria curvas cota-volume-vazão, rendimento,
defluências, cascatas e restrições ambientais. O modelo reduz esses elementos a um
reservatório equivalente por subsistema. Como $h=P^h Delta t$ e a afluência $A$
já está convertida em MWh-equivalentes, a dinâmica é

$ S_(d,ell,t+1)=S_(d,ell,t)+A_(d,ell,t)-h_(d,ell,t)-v_(d,ell,t). $

As bandas e o fechamento são

$ 0 <= S_(d,ell,t) <= overline(S)_ell, quad
  0 <= h_(d,ell,t) <= overline(H)_ell, quad
  v_(d,ell,t) >= 0, $

$ S_(d,ell,24)=S_(d,ell,0)=0,5 overline(S)_ell. $

A capacidade energética equivalente é definida como oito horas da potência hídrica
proxy. A condição cíclica impede que um medoide gere energia esvaziando um estoque
inicial gratuito. Ela permite arbitragem intradiária, mas não cria uma ligação
sazonal fictícia entre medoides.

#imgfig(
  [Estado do reservatório, afluência, turbinamento e vertimento.],
  "outputs/figures/rel_reservoir_state.png",
  note: [
    *Como ler.* O painel superior acompanha $S_t$ do instante 0 ao 24 em um medoide
    do Sudeste/CO. O painel inferior coloca na mesma escala os três termos de fluxo:
    $A_t$, $h_t$ e $v_t$. *Resultado.* A curva retorna exatamente ao estoque
    inicial e o erro máximo de $S_(t+1)-S_t-A_t+h_t+v_t$ é
    $5,47 times 10^(-11)$ MWh. Isso é uma auditoria da restrição, não apenas uma
    ilustração. *Cuidado.* Estoque em MWh-equivalentes não deve ser interpretado
    como volume físico de um reservatório específico.
  ]
)

== Problema de despacho

Se $p_d$ é o número de dias representado pelo medoide, o termo principal da função
objetivo é

$
  min sum_d p_d sum_(ell,t) [
    sum_b (c_(ell,b)^"th" n_(d,ell,b,t)
      + frac(epsilon_n,2) n_(d,ell,b,t)^2)
    + pi_h h_(d,ell,t)
    + pi_s c_(d,ell,t)^s
    + pi_w c_(d,ell,t)^w
    + pi_u u_(d,ell,t)
    + gamma_r (n_(d,ell,t)-n_(d,ell,t-1))^2
  ] + epsilon_x ||x||_2^2.
$

O vetor $x$ reúne as variáveis de despacho regularizadas. No baseline,
$pi_u=8.327,76$ R\$/MWh, $pi_h=250$ R\$/MWh, $pi_s=pi_w=0,01$
R\$/MWh, $epsilon_n=2 times 10^(-4)$, $gamma_r=2 times 10^(-3)$ e
$epsilon_x=10^(-2)$. As penalidades de corte são epsilons de desempate, não uma
representação de compensação regulatória.

O balanço nodal usa a convenção

$
  L_(d,ell,t)
  - g_(d,ell,t)^s - g_(d,ell,t)^w - h_(d,ell,t)
  - sum_b n_(d,ell,b,t) - g_(d,ell,t)^"nuc"
  - F_(d,ell,t)^"net" - u_(d,ell,t) = 0.
$

Os blocos térmicos satisfazem $0<=n_(ell,b,t)<=overline(N)_(ell,b)$, os
fluxos satisfazem $-overline(F)_e<=F_(e,t)<=overline(F)_e$ e o hub satisfaz
soma algébrica nula.

== Dual, teorema do envelope e KKT

Se $tilde(lambda)_(d,ell,t)$ é o dual bruto do balanço e o objetivo do medoide é
ponderado por $p_d$, o preço reportado é normalizado como

$ lambda_(d,ell,t)=frac(tilde(lambda)_(d,ell,t),p_d). $

Com o balanço escrito como carga menos geração igual a zero, o teorema do envelope
fornece diretamente

$ lambda_(d,ell,t)=frac(partial V^*, partial L_(d,ell,t)). $

Não há troca posterior de sinal. Em um bloco térmico interior e omitindo apenas
os termos de banda, a condição de primeira ordem é

$ lambda_t = c_b^"th" + epsilon_n n_(b,t)
  + 2 gamma_r(2n_t-n_(t-1)-n_(t+1))
  + "termo do regularizador". $

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

Mantendo a afluência de referência, o notebook varia $pi_h$ entre 80, 250 e 500
R\$/MWh. A média de $lambda$ é ponderada pelos pesos dos medoides, assim como os
totais anuais.

#text(size: 8.7pt)[
#table(
  columns: (1.02fr, 1.05fr, 1.12fr, 1.12fr, 1.12fr),
  inset: 4pt,
  align: right,
  [$pi_h$ R\$/MWh], [$lambda$ médio], [Hidro TWh], [Térmica TWh], [Vertimento TWh],
  [80], [275,83], [385,57], [85,19], [10,02],
  [250], [321,35], [358,00], [112,76], [37,59],
  [500], [465,16], [335,27], [135,50], [60,33],
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
  [Seco -30%], [5.834,95], [134,59], [264,48], [12,43], [71,69],
  [Base], [321,35], [112,76], [358,00], [37,59], [0,00],
  [Úmido +30%], [241,00], [86,97], [383,80], [130,47], [0,00],
)
]

#imgfig(
  [Sensibilidade hidrológica com as demais entradas fixas.],
  "outputs/figures/rel_hydrology_sensitivity.png",
  note: [
    *Como ler.* O painel esquerdo mostra o preço-sombra médio anual ponderado; o
    centro mostra déficit anualizado; o direito separa térmica, hidro e vertimento.
    *Resultado.* Reduzir afluência em 30% cria 71,69 TWh de déficit e eleva o preço
    médio a 5.834,95 R\$/MWh, ainda abaixo do VOLL porque é uma média sobre nós e
    horas. No cenário úmido, a térmica cai para 86,97 TWh e o vertimento sobe para
    130,47 TWh. *Cuidado.* O choque é uniforme e determinístico. Ele não representa
    uma distribuição de ENA nem dependência entre bacias, e o vertimento é do
    reservatório equivalente.
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

Antes de tomar o limite de muitos agentes, isolamos o mecanismo com dois projetos
solares no Sudeste/CO. O produtor 1 é o incumbente, com
$K_1=11.446,23$ MW. O produtor 2 entra com

$ K_2=kappa K_1, quad kappa in (0;0,25;dots;1,50). $

Ambos enfrentam o mesmo clearing com rede, térmica, hidro, reservatório, eólica,
nuclear e déficit. Se a geração potencial de cada produtor é

$ overline(g)_(i,d,t)^s=K_i xi_i a_(d,"SE",t)^s, $

a fração efetivamente entregue pelo nó é aplicada proporcionalmente aos dois
projetos quando há corte agregado. A receita anual e o preço capturado são

$ Rev_i=sum_d p_d sum_t lambda_(d,"SE",t) g_(i,d,t)^s, $

$ P_i^"cap"=frac(Rev_i, sum_d p_d sum_t g_(i,d,t)^s). $

O lucro usa o custo anualizado do tipo:

$ Pi_i=Rev_i-q_("SE",theta_i)K_i. $

Quando $kappa=0$, o entrante não existe: sua receita é zero, mas seu preço
capturado é indefinido porque o denominador também é zero. A tabela e a figura não
o substituem por um zero artificial.

#imgfig(
  [Incumbente e entrante: preço capturado, lucro e curtailment.],
  "outputs/figures/rel_finite_two_player.png",
  note: [
    *Como ler.* O eixo horizontal é $kappa=K_2/K_1$. O primeiro painel compara
    preço capturado; o segundo, lucro anual; o terceiro, corte alocado a cada
    projeto. *Resultado.* Entre $kappa=0$ e $0,50$, o preço capturado do incumbente
    cai de 3.011,23 para 2.967,37 R\$/MWh enquanto o déficit anualizado cai de 1,615
    TWh para 0,00183 TWh. Entre $kappa=0,50$ e $0,75$, o sistema cruza um limiar de
    adequação: o déficit desaparece e o preço capturado cai abruptamente para
    369,03 R\$/MWh. Em $kappa=1,50$, ele chega a 313,23 R\$/MWh. O curtailment
    solar permanece numericamente desprezível nessa faixa. *Cuidado.* A queda
    descontínua é um efeito de mudança de restrição ativa, especialmente saída do
    déficit, e não prova uma função de preço suave ou uma estratégia de Nash entre
    dois agentes grandes.
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

Para uma localização $ell$ e tipo $theta$, a receita spot anual por MW gerada pelo
clearing é

$ R_(ell,theta)(m)=sum_d p_d sum_t
  lambda_(d,ell,t)(m) a_(d,ell,t)^s xi_theta
  delta_(d,ell,t)^"delivery"(m), $

onde $delta^"delivery"=g^s/overline(g)^s$ incorpora curtailment. A dependência em
$m$ é essencial: a distribuição determina capacidade, a capacidade determina o
clearing e o clearing determina a receita.

#imgfig(
  [Receita anual por MW como função da capacidade solar agregada.],
  "outputs/figures/canibalizacao_R_de_C.png",
  note: [
    *Como ler.* O eixo horizontal multiplica simultaneamente a capacidade solar de
    referência nas quatro regiões; cada curva é a receita anual por MW de um
    projeto naquela localização, calculada novamente pelo clearing. *Resultado.*
    As curvas caem com capacidade, confirmando no operador completo a
    realimentação de canibalização derivada na seção empírica. As inclinações e os
    níveis diferem por rede, perfil renovável e disponibilidade de recursos.
    *Cuidado.* Monotonicidade observada nessa grade favorece estabilidade, mas não
    demonstra que o operador de melhor resposta seja uma contração global.
  ]
)

== Benchmark estático de campo médio

O benchmark reduzido mantém a massa regional fixa e representa uma escolha
contínua de capacidade. Para cada $(ell,theta)$,

$ U_(ell,theta)(k)
  = (R_(ell,theta)-q_(ell,theta)+p) kappa_u k
  - frac(gamma,2) k^2, $

onde $p$ é um prêmio não spot opcional, $kappa_u$ converte o índice de capacidade
em MW e $gamma$ controla curvatura. Com $tau=0$, o ótimo contínuo é

$ k^* = pos(frac((R-q+p)kappa_u,gamma)), $

truncado ao intervalo da grade. Quando $k^*$ fica entre os nós $b$ e $b+1$, sua
representação baricêntrica usa pesos $1-phi$ e $phi$, com
$phi=k^*-b$. Essa representação preserva exatamente o momento de capacidade;
ela não é chamada de mistura ótima de um jogo discreto, pois bins adjacentes não
precisam ser indiferentes.

O ponto fixo é

$ m=BR(R(m)). $

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
das respectivas tolerâncias.

== Calibração por entrada livre

O nível agregado é calibrado para a proxy observada de 39,069 GW. Com a curvatura
inicial e prêmio zero, o modelo escolhe 41,37 GW; por isso a bisseção aumenta
$gamma$ até 0,443 milhão de R\$ por bin ao quadrado e mantém prêmio não spot igual
a zero. A solução final contém 39,088 GW.

#warnbox([O que a calibração não prova])[ 
  Reproduzir o total observado é uma condição de calibração, não uma validação fora
  da amostra. A distribuição regional, a sensibilidade a custos, os preços e a
  dinâmica continuam sendo resultados condicionais. Como a massa regional fica
  fixa nessa camada, ela escolhe capacidade dentro da localização atribuída; a
  migração endógena pertence ao MFG dinâmico.
]

#imgfig(
  [Convergência e capacidade do benchmark estático.],
  "outputs/figures/rel_mfg_convergencia.png",
  note: [
    *Como ler.* À esquerda, escalas logarítmicas comparam $r_C$, $r_m$ e o passo
    relaxado. À direita, a capacidade agregada converge para o alvo de calibração.
    *Resultado.* A rotina termina na iteração 146 com
    $r_C=5,25 times 10^(-5)$, $r_m=4,75 times 10^(-4)$ e
    $r_R=3,65 times 10^(-7)$. Todos satisfazem as tolerâncias e o ponto é
    certificado. *Cuidado.* A curva do passo fica abaixo dos resíduos por um fator
    ligado a $alpha$; somente os resíduos do operador sustentam o rótulo de ponto
    fixo.
  ]
)

#text(size: 8.7pt)[
#table(
  columns: (0.72fr, 1.00fr, 1.00fr, 1.10fr),
  inset: 4pt,
  align: right,
  [Subs.], [Proxy observada GW], [Privado GW], [Privado - proxy GW],
  [SE], [19,077], [24,972], [+5,895],
  [S], [5,919], [4,403], [-1,516],
  [NE], [11,789], [7,642], [-4,148],
  [N], [2,284], [2,072], [-0,213],
  [*Total*], [*39,069*], [*39,088*], [*+0,019*],
)
]

A calibração fecha o total, mas o campo de receitas redistribui capacidade dentro
das massas regionais: o Sudeste/CO cresce, enquanto Sul, Nordeste e Norte ficam
abaixo das proxies. Isso é uma implicação interna do modelo, não uma recomendação
de localização.

== Contrafactual de total fixo

Para preservar a ideia de “orçamento fixo” sem chamá-la de equilíbrio privado,
resolvemos o problema do planejador com a restrição

$ sum_ell K_ell=K^"obs"=39,069 " GW". $

Capacidade total e todas as demais entradas ficam fixas; apenas a localização é
otimizada junto com o despacho. O resultado concentra 34,753 GW no Sudeste/CO e
4,316 GW no Nordeste, com valores numericamente nulos no Sul e no Norte. Essa
solução de canto mostra que, sob os custos e a rede reduzida, o operador valoriza
fortemente o Sudeste/CO. Ela deve ser lida como diagnóstico de identificação: uma
alocação tão concentrada sinaliza que custos locacionais, limites de rede e
potenciais máximos ainda precisam de calibração externa antes de uso normativo.

== Planejador social com total livre e limites regionais

O planejador escolhe capacidade e despacho simultaneamente:

$
  min_(K,x) sum_ell q_ell K_ell + C^"op"(K,x)
$

sujeito ao mesmo clearing e às bandas

$ 0<=K_ell<=3K_ell^"obs". $

O benefício marginal da capacidade, quando interior, decorre do envelope:

$ frac(partial C^"op",partial K_ell)
  =-sum_d p_d sum_t lambda_(d,ell,t) a_(d,ell,t)^s
    delta_(d,ell,t)^"delivery" + "termos de corte". $

Logo, a condição interior iguala custo marginal de investimento e redução
marginal do custo sistêmico. O agente privado toma $R(m)$ como dado; o planejador
internaliza o efeito agregado sobre déficit, congestionamento, água, térmica e
curtailment.

#imgfig(
  [Capacidade observada, total fixo, privada e social por subsistema.],
  "outputs/figures/rel_mercado_vs_social.png",
  note: [
    *Como ler.* Cinza é a proxy `p99`; laranja é o contrafactual com total fixo;
    verde é o benchmark privado calibrado; azul é o planejador com total livre e
    teto regional. *Resultado.* Os três primeiros totais são aproximadamente 39,1
    GW por construção ou calibração. O planejador escolhe 77,91 GW: 57,23 no
    Sudeste/CO, 12,31 no Nordeste, 5,57 no Sul e 2,80 no Norte. *Restrição ativa.*
    O valor do Sudeste/CO é exatamente três vezes sua proxy, portanto o teto
    regional está ativo; 77,91 GW é o ótimo do problema limitado, não um ótimo
    irrestrito em capacidade. *Interpretação.* Nesta parametrização, o custo social
    evitado, sobretudo por adequação, supera o custo de nova solar até níveis acima
    do benchmark privado. Isso inverte o resultado preliminar do relatório antigo
    e deve ser tratado como resultado de modelo, não como recomendação de dobrar a
    capacidade.
  ]
)

#text(size: 8.35pt)[
#table(
  columns: (0.70fr, 0.95fr, 1.10fr, 0.95fr, 1.05fr),
  inset: 3.8pt,
  align: right,
  [Subs.], [Proxy GW], [Total fixo GW], [Privado GW], [Social limitado GW],
  [SE], [19,077], [34,753], [24,972], [57,231],
  [S], [5,919], [0,000], [4,403], [5,568],
  [NE], [11,789], [4,316], [7,642], [12,306],
  [N], [2,284], [0,000], [2,072], [2,801],
  [*Total*], [*39,069*], [*39,069*], [*39,088*], [*77,906*],
)
]

#resultbox([Leitura correta da comparação])[ 
  A comparação não identifica empiricamente “subinvestimento” ou “sobreinvestimento”.
  Ela mostra que, dentro do modelo calibrado, o agente privado e o planejador
  valorizam margens diferentes. O resultado social é sensível ao VOLL, ao custo
  anualizado, ao teto regional e à representação de contratos. O fato de um teto
  estar ativo deve motivar sensibilidade adicional antes de qualquer conclusão
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

Em um jogo finito com muitos projetos, o agente $i$ resolveria

$ max_(a_i) J_i(a_i,a_(-i)). $

O campo médio substitui o vetor enorme $a_(-i)$ pela distribuição $m_y$. Cada
agente é infinitesimal e toma $(m,R)$ como dado ao resolver seu Bellman; a
população agregada, contudo, determina o mesmo $(m,R)$ que os agentes antecipam.

== Estado, ação e transição

O estado é

$ x=(ell,theta,b) in cal(X)=cal(L) times cal(Theta) times cal(B), $

onde $b in (0,dots,12)$ e cada passo corresponde a 93 MW por agente. A ação é

$ a=(ell',delta), quad delta in (-1,0,1), $

e a transição determinística é

$ T((ell,theta,b),(ell',delta))
  =(ell',theta,b'), quad b'=b+delta, $

para ações que mantêm $b'$ entre 0 e 12. O tipo $theta$ é permanente. A
distribuição satisfaz

$ m_y(ell,theta,b)>=0, quad sum_(ell,theta,b)m_y(ell,theta,b)=1. $

Se $M$ é o número equivalente de agentes, a capacidade no nó é

$ K_(ell,y)(m_y)=M sum_(theta,b) K_b m_y(ell,theta,b),
  quad K_b=93b " MW". $

O número $M$ é escolhido para que o estado inicial, concentrado em $b=4$, reproduza
os 39,069 GW da proxy observada. Isso fixa a escala da população, não o caminho
posterior.

== Recompensa individual

Dada a trajetória de receitas anuais por MW produzida pelo clearing, a recompensa
de escolher $(ell',b')$ a partir de $(ell,b)$ é

$
  r_y(x,a;m_y,R_y)
  = (R_(y,ell',theta)+p-q_(ell',theta))K_(b')
  - frac(chi,2)(K_(b')-K_b)^2
  - F^"entry" ind_(b=0 " e " b'>0)
  - rho^"move" ind_(ell' != ell)
  - gamma^"conc" s_(ell')(m_y)(K_(b')/1000)^2,
$

onde

$ s_ell(m_y)=sum_(theta,b)m_y(ell,theta,b). $

Os parâmetros dinâmicos são $H=6$, $beta=0,95$, $chi=2.000$
R\$/MW², $F^"entry"=20$ milhões de R\$, $rho^"move"=500$ milhões de R\$ e
$gamma^"conc"=5$ milhões de R\$. O prêmio $p$ vem da calibração estática e é
zero no baseline.

O custo de concentração é reduzido e não deve ser interpretado como restrição de
planejamento territorial observada. Ele evita que pequenas diferenças de receita
movam toda a massa para um único nó no modelo agregado.

== Bellman regularizado e política

Defina

$ Q_y(x,a)=r_y(x,a;m_y,R_y)+beta V_(y+1)(T(x,a)). $

Com $tau>0$, o valor regularizado é

$ V_y(x)=tau log sum_(a in cal(A)(x)) exp(Q_y(x,a)/tau), $

e a política é

$ pi_y(a|x)=frac(exp(Q_y(x,a)/tau),
  sum_(a' in cal(A)(x)) exp(Q_y(x,a')/tau)). $

O baseline usa $tau=100$ milhões de R\$. Portanto, o objeto calculado é
explicitamente um *MFG entrópico regularizado*. A aleatoriedade não é acrescentada
depois da otimização: ela é a política ótima do problema regularizado.

Quando $tau=0$, o código substitui o `log-sum-exp` pelo máximo de Bellman e divide
probabilidade apenas entre ações numericamente empatadas:

$ pi_y(a|x)>0 arrow.r Q_y(x,a)=max_(a')Q_y(x,a'). $

Não se afirma que o equilíbrio regularizado com $tau>0$ seja idêntico ao jogo
discreto não regularizado.

== Equação forward

Dada a política, a distribuição evolui por

$
  hat(m)_(y+1)(x')
  =sum_(x in cal(X))sum_(a in cal(A)(x))
    m_y(x) pi_y(a|x) ind_(T(x,a)=x').
$

Como a transição é determinística e cada política soma um, a equação conserva
massa. Em memória, o notebook verifica numericamente

$ max_y abs(sum_x m_y(x)-1)=2,22 times 10^(-16), $

e, para cada estado e estágio,

$ max_(y,x) abs(sum_a pi_y(a|x)-1)=4,44 times 10^(-16). $

Depois da serialização em CSV e nova agregação, os limites conservadores observados
são $4,44 times 10^(-15)$ para massa e $1,11 times 10^(-15)$ para política.

== Clearing e condição de ponto fixo

Em cada estágio $y$, $m_y$ produz $K_(ell,y)$; o clearing compartilhado devolve
$R_(y,ell,theta)$; o Bellman devolve $pi_y$; e a equação forward devolve
$hat(m)$. Em notação de operador,

$ m=cal(F)(BR(cal(R)(m))), $

onde a igualdade vale para a trajetória $m=(m_0,dots,m_H)$ inteira, com $m_0$
fixo.

A iteração usa

$ m^(n+1)=(1-alpha)m^n+alpha hat(m)^n, quad alpha=0,10, $

mas a certificação é feita no operador sem relaxação:

$ r_m=||hat(m)-m||_infinity, $

$ r_C=max_(y,ell) frac(abs(hat(K)_(y,ell)-K_(y,ell)),K_(y,ell)+1), $

$ r_R=max_(y,ell,theta)
  frac(abs(hat(R)_(y,ell,theta)-R_(y,ell,theta)),
       abs(R_(y,ell,theta))+1). $

As tolerâncias são $2 times 10^(-4)$ para $r_m$ e
$5 times 10^(-3)$ para $r_C$ e $r_R$.

#imgfig(
  [Certificação do MFG dinâmico e trajetória de capacidade.],
  "outputs/figures/rel_dynamic_mfg_fixed_point.png",
  note: [
    *Como ler.* O painel esquerdo usa escala logarítmica para os resíduos do
    operador antes da relaxação e para o passo numérico. O painel direito mostra
    $K_(ell,y)$ em cada estágio. *Resultado numérico.* A iteração 67 termina com
    $r_m=1,91 times 10^(-4)$, $r_C=3,90 times 10^(-4)$ e
    $r_R=5,85 times 10^(-6)$; os três passam suas tolerâncias. *Resultado
    econômico.* A capacidade total parte de 39,069 GW e chega a 55,376 GW, com
    expansão concentrada no Sudeste/CO. *Cuidado.* Convergência comprova
    consistência interna do MFG regularizado, não acurácia empírica dos parâmetros.
  ]
)

== Trajetória agregada

#text(size: 8.35pt)[
#table(
  columns: (0.72fr, 1.02fr, 1.02fr, 1.02fr, 1.02fr),
  inset: 3.8pt,
  align: right,
  [Subs.], [Massa inicial], [Massa terminal], [Cap. inicial GW], [Cap. terminal GW],
  [SE], [48,83%], [60,08%], [19,077], [35,729],
  [S], [15,15%], [11,08%], [5,919], [5,278],
  [NE], [30,18%], [22,16%], [11,789], [10,978],
  [N], [5,85%], [6,67%], [2,284], [3,392],
  [*Total*], [*100%*], [*100%*], [*39,069*], [*55,376*],
)
]

A capacidade pode aumentar mesmo quando a massa de uma região cai, pois os agentes
remanescentes podem migrar para bins maiores. Esse é o motivo para reportar massa e
capacidade separadamente.

#imgfig(
  [Trajetória completa da distribuição por localização e capacidade individual.],
  "outputs/figures/rel_dynamic_mfg_distribution_trajectory.png",
  note: [
    *Como ler.* Cada painel corresponde a um subsistema. O eixo horizontal é o
    estágio $y=0,dots,6$; o vertical é a capacidade individual $K_b$; a cor é
    $sum_theta m_y(ell,theta,b)$, isto é, massa incondicional da população naquele
    nó e bin. A escala de cor é comum aos quatro painéis. *Resultado.* A massa
    inicial está concentrada em 372 MW por agente. Ao longo do horizonte, o
    Sudeste/CO desloca massa para bins de 558, 651, 744 MW e superiores; Nordeste,
    Sul e Norte também se espalham, mas com menor peso agregado. *Cuidado.* Uma
    célula escura pode significar pouca massa regional ou dispersão por muitos
    bins. Por isso a figura deve ser lida junto com a massa locacional abaixo.
  ]
)

#imgfig(
  [Evolução da fração da população por subsistema.],
  "outputs/figures/rel_dynamic_mfg_location_mass.png",
  note: [
    *Como ler.* A altura total da área empilhada é sempre um; cada cor integra
    tipos e bins de capacidade dentro de uma localização. *Resultado.* A parcela
    do Sudeste/CO cresce de 48,83% para 60,08%, enquanto Nordeste e Sul perdem
    massa relativa. O Norte cresce levemente. *Auditoria.* A soma das áreas em cada
    estágio difere de um por no máximo $2,22 times 10^(-16)$. *Cuidado.* A figura
    mostra realocação de agentes, não fluxo físico de energia entre subsistemas.
  ]
)

== Política de localização e capacidade

#imgfig(
  [Política no estágio inicial: realocação e ajuste esperado.],
  "outputs/figures/rel_dynamic_mfg_policy.png",
  note: [
    *Como ler.* As duas colunas separam tipos `high` e `low`. Na linha superior, a
    cor é a probabilidade total de escolher $ell'!=ell$; na inferior, é
    $EE[Delta K|x]$ em MW, com vermelho para expansão e azul para redução. As
    linhas são localizações atuais e as colunas são bins atuais. *Resultado.*
    Estados do Sudeste/CO têm probabilidade de realocação próxima de zero, mas
    ajuste esperado positivo em quase toda a grade. Estados do Nordeste e do Sul
    apresentam maior propensão a migrar nos bins altos. No bin máximo, o ajuste
    esperado se torna negativo porque a ação $+1$ é inviável e a política ainda
    atribui massa a redução. *Cuidado.* A regularização entrópica suaviza escolhas;
    cores intermediárias não são erro numérico nem heterogeneidade observada.
  ]
)

#imgfig(
  [Trajetórias de ações modais para agentes representativos.],
  "outputs/figures/rel_dynamic_mfg_modal_paths.png",
  note: [
    *Como ler.* Cada linha fixa a localização inicial, cada coluna é um estágio e
    cada célula informa o destino modal e a capacidade após a ação modal. Os dois
    painéis separam tipos. *Uso correto.* A figura oferece uma narrativa individual
    intuitiva para a política. *Limite fundamental.* A ação modal descarta todas as
    probabilidades não modais; somar essas trajetórias não reconstrói $m_y$. Para
    qualquer conclusão agregada, valem a equação forward e os mapas de distribuição,
    não a trajetória modal.
  ]
)

#resultbox([Resultado central do MFG dinâmico])[ 
  Existe uma trajetória internamente consistente do MFG entrópico regularizado nas
  tolerâncias declaradas. Ela desloca massa para o Sudeste/CO, aumenta o bin médio
  de capacidade e leva o total de 39,069 a 55,376 GW. O resultado é diferente do
  benchmark estático porque agora localização, capacidade e valor futuro são
  escolhidos em conjunto. O rótulo “certificado” refere-se ao ponto fixo do modelo,
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
  [Reservatório], [$5,47 times 10^(-11)$ MWh], [erro máximo da equação de movimento; fechamento cíclico exato.],
  [Cenário seco], [71,69 TWh de déficit], [choque uniforme de -30% na afluência, demais entradas fixas.],
  [Benchmark estático], [39,088 GW], [calibrado para 39,069 GW e certificado por três resíduos.],
  [Total fixo], [39,069 GW], [contrafactual de realocação; não é equilíbrio privado.],
  [Social limitado], [77,906 GW], [teto do Sudeste/CO ativo; resultado condicional, não irrestrito.],
  [MFG dinâmico], [55,376 GW terminal], [MFG entrópico regularizado certificado na iteração 67.],
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
- A rede tem cinco nós contábeis e quatro interfaces fixas; perdas, contingências e
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

- O custo locacional e a produtividade dos dois tipos são proxies. O social atingir
  o teto regional mostra que a fronteira de capacidade importa materialmente.
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

A versão 2 preserva a pergunta, a arquitetura e a matemática central do relatório
original, mas transforma aproximações preliminares em objetos verificáveis. Carga
residual e captura explicam o mecanismo; o clearing reúne restrições físicas em um
único primal; KKT dá interpretação ao preço-sombra e ao valor da água; o laboratório
finito mostra como entrada altera restrições ativas; e o MFG fecha Bellman,
distribuição forward e clearing em uma trajetória consistente.

As principais correções são substantivas: curtailment explícito, estoque hídrico
com fechamento cíclico, sinal do dual determinado pela convenção do balanço,
distinção entre passo relaxado e resíduo do operador, representação baricêntrica sem
falsa mistura discreta, planejador independente e nomenclatura explícita do MFG
entrópico. Os resultados atuais não precisam coincidir com os antigos porque os
objetos computacionais mudaram; o cerne econômico e matemático, entretanto, é o
mesmo e agora está mais rigoroso.

#resultbox([Mensagem final])[ 
  A conclusão defensável não é um número isolado de capacidade. É a existência de
  uma cadeia auditável: dados horários determinam perfis; perfis e capacidade
  determinam o clearing; o clearing determina receitas e preços-sombra; agentes
  respondem a esses sinais; e a distribuição resultante deve reproduzir o campo
  usado para calculá-los. Quando essa última igualdade é afirmada, ela vem
  acompanhada do resíduo que a certifica.
]

#pagebreak(weak: true)

= Apêndice: rastreabilidade computacional

== Notebooks

#table(
  columns: (1.25fr, 3.75fr),
  inset: 4pt,
  align: left,
  [`06_implementacao_documento_mfg_v4.ipynb`], [referência histórica consultada; mantém a versão anterior do clearing e da calibração.],
  [`07_mfg_resumo.ipynb`], [fonte do modelo, parâmetros, execução, tabelas, auditorias e certificados.],
  [`08_relatorio_figuras.ipynb`], [camada de visualização; consome o painel e as tabelas do 07 sem reestimar o modelo.],
)

== Tabelas de auditoria

#text(size: 8.55pt)[
#table(
  columns: (1.85fr, 3.15fr),
  inset: 4pt,
  align: left,
  [`reservoir_state_path.csv`], [trajetórias horárias de estoque, afluência, hidro e vertimento.],
  [`reservoir_audit.csv`], [erro da dinâmica e diferença entre estoque inicial e final.],
  [`hydrology_sensitivity.csv`], [seco, base e úmido com as demais entradas fixas.],
  [`finite_two_player_sensitivity.csv`], [receita, captura, lucro, corte, preço e déficit por $kappa$.],
  [`mfg_convergence_history.csv`], [passos e resíduos do benchmark estático.],
  [`mfg_locational_equilibrium.csv`], [capacidade, receita, custo e margem por localização e tipo.],
  [`fixed_budget_capacity.csv`], [realocação com capacidade total fixada na proxy.],
  [`private_vs_social_capacity.csv`], [proxy, total fixo, privado e social por subsistema.],
  [`dynamic_mfg_convergence.csv`], [histórico de $r_m$, $r_C$ e passo relaxado.],
  [`dynamic_mfg_distribution.csv`], [trajetória longa de $m_y(ell,theta,b)$.],
  [`dynamic_mfg_policy.csv`], [uma linha por ação com $pi_y(a|x)$.],
  [`dynamic_mfg_policy_summary.csv`], [momentos e ação modal por estado.],
  [`dynamic_mfg_path.csv`], [massa, capacidade, receita e clearing por estágio e local.],
  [`dynamic_mfg_modal_paths.csv`], [trajetórias modais auxiliares para leitura microeconômica.],
)
]

== Checklist de reprodução

#enum[
  executar `00_construcao_painel.ipynb` quando os caches de dados forem refeitos;
][
  executar integralmente `07_mfg_resumo.ipynb` e exigir zero célula com erro;
][
  verificar os `asserts` de balanço hídrico, fechamento cíclico, massa e política;
][
  conferir `certified=True` e os resíduos não relaxados nas duas rotinas;
][
  executar `08_relatorio_figuras.ipynb` para regenerar as imagens do relatório;
][
  compilar este arquivo e revisar visualmente todas as páginas, legendas e tabelas.
]

#v(1em)
#align(center)[
  #text(size: 8.3pt, fill: gray.darken(20%))[
    Versão 2 gerada a partir dos notebooks 07 e 08 executados integralmente.
  ]
]
