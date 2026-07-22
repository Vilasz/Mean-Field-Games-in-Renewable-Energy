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


= Dados, notação e compressão temporal

== Índices, estados e unidades

O tempo operacional é horário, $t=0,dots,23$. Os dias representativos são
$d in cal(D)$, os subsistemas são

$ cal(L) = ("SE", "S", "NE", "N"), $

os blocos térmicos são $b in cal(B)^"th"$ e os tipos solares são
$theta in ("high", "low")$. No MFG dinâmico, $y=0,dots,H$ denota o estágio de
decisão e não a hora do despacho.

*O que são os blocos térmicos.* A geração térmica do SIN reúne dezenas de usinas
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

*Por que os blocos importam.* Eles são o mecanismo por trás de três resultados
centrais do relatório. Primeiro, o preço-sombra: nas horas em que a térmica é a
tecnologia marginal, o bloco acionado no topo da pilha fixa $lambda$, de modo que a
condição de primeira ordem $lambda_t = c_b^"th" + dots$ da seção do clearing liga
diretamente preço a mérito térmico. Segundo, a adequação: solar e eólica não
produzem nas horas sem sol nem sob rampa do entardecer, e são os blocos térmicos
(com a hidro) que preenchem essa carga residual; sem uma pilha térmica com potência
suficiente, o problema recorre ao déficit penalizado. Terceiro, a canibalização: ao
crescer, a solar desloca os blocos térmicos *mais baratos* nas horas de sol,
empurrando a tecnologia marginal — e o preço — para baixo justamente quando os
projetos solares vendem. Sem uma representação explícita do mérito térmico, esse
canal de preço teria de ser imposto por uma curva arbitrária.

*O que são os tipos solares.* Nem todo projeto solar é igual: sítios diferem em
irradiação, fator de capacidade e custo de instalação. Em vez de tratar cada
projeto individualmente, resumimos essa heterogeneidade em dois *tipos*
$theta in ("high", "low")$. O tipo `high` tem produtividade relativa
$xi_"high"=1,08$ com multiplicador de custo 0,94 (mais energia por MW, mais barato);
o tipo `low` tem $xi_"low"=0,94$ com multiplicador 1,08 (menos energia por MW, mais
caro). As massas iniciais são 45% `high` e 55% `low`. O fator $xi_theta$ entra
diretamente na geração disponível por MW instalado, de modo que um mesmo MW rende
receitas diferentes conforme o tipo.

*Por que os tipos importam.* No benchmark estático, os tipos separam a decisão de
investimento por qualidade de recurso: a entrada livre e o custo anualizado
$q_(ell,theta)$ dependem de $theta$, e a distribuição por tipo determina quanto da
capacidade é rentável em cada localização. No MFG dinâmico, o tipo é uma *dimensão
permanente do estado* $x=(ell,theta,b)$ — um agente não muda de tipo ao migrar ou
ao mudar de bin —, o que torna a distribuição populacional $m_y$ e a condição de
consistência de campo médio genuinamente heterogêneas. Reduzir tudo a um tipo médio
apagaria essa seleção e daria uma resposta de investimento artificialmente uniforme.

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

*O problema: compressão temporal.* Cada camada de decisão do relatório resolve o
clearing muitas vezes — o MFG dinâmico, por exemplo, o resolve a cada estágio e a
cada iteração do ponto fixo. Rodar os 365 despachos de 24 horas do ano em cada uma
dessas chamadas seria desnecessariamente caro sem alterar as conclusões, pois muitos
dias têm perfis horários de carga e de geração renovável quase idênticos. A ideia é
substituir o calendário completo por um punhado de *dias representativos*, cada um
carregando um peso igual ao número de dias reais que ele resume.

*O que é um cluster de dias representativos.* Descrevemos cada dia do ano por um
vetor de características (perfil horário de carga, de solar e de eólica, além de
marcadores de estação e de calendário) e agrupamos os 365 vetores por similaridade.
Cada grupo — um *cluster* — reúne dias com comportamento operacional parecido: por
exemplo, dias úteis de inverno com pouca solar, ou fins de semana de verão com forte
geração diurna. Em vez de simular todos os dias do cluster, simulamos apenas um dia
que o represente e anexamos a ele o peso do grupo. O número de clusters é escolhido
entre 4 e 12 combinando três critérios: a *silhueta* (uma medida de quão bem
separados e coesos ficam os grupos), a cobertura das quatro *estações* e a distinção
entre *regimes de calendário* (útil, fim de semana, feriado). A tabela abaixo usa
quatro blocos.

*O que é o medoide e por que o escolhemos.* O representante de cada cluster é o seu
*medoide*: o dia real cujo perfil está mais próximo dos demais dias do grupo, isto
é, o que minimiza a soma das distâncias aos outros membros. Ele contrasta com o
*centroide* do $k$-médias, que é a média ponto a ponto do grupo — um perfil
*sintético* que pode não corresponder a nenhum dia observado e que suaviza
artificialmente rampas e picos ao promediar dias ligeiramente defasados. Como o
clearing depende criticamente da forma horária (o vale solar do meio-dia, a rampa do
entardecer, o horário do pico de carga), preferimos o medoide por três razões: ele
preserva um perfil horário genuinamente observado, sem inventar transições que nunca
ocorreram; é robusto a dias atípicos, que não distorcem um representante real como
distorceriam uma média; e mantém a interpretação física — cada bloco tem uma data
real associada, listada na tabela.

*As categorias de calendário e por que elas entram.* Além da estação, cada dia
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
    operacionais já saem anualizados, sem multiplicar de novo por 365. O que a figura
    também revela, e convém não esconder, é o custo dessa compressão: quatro blocos
    são poucos para o Brasil. Um único dia de verão carrega 71 dias reais e um único
    dia de inverno carrega 104, então toda a variabilidade dentro de cada estação —
    ondas de calor, frentes frias, dias nublados atípicos — colapsa em um perfil só.
    Além disso, os medoides são dias isolados: eles preservam a forma de cada dia,
    mas não preservam a sequência entre dias, isto é, não há memória de uma seca que
    se agrava por semanas nem de um reservatório enchendo. É por isso que, mais
    adiante, cada reservatório é fechado ciclicamente dentro do bloco, uma escolha
    que evita ganhar energia de graça mas que também impede representar transferência
    sazonal de água — uma limitação real, não um detalhe técnico.
  ]
)

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
    que mais chama atenção é quão *fraca* é essa relação: mesmo no melhor caso, a
    carga residual explica pouco da variação horária do preço, sinal de que água,
    déficit e rede pesam tanto quanto o nível residual. O Nordeste é o caso
    instrutivo: a correlação é negativa, -0,122, porque a eólica local muitas vezes
    excede a carga regional e o preço passa a ser ditado por excedente e intercâmbio,
    não pela carga própria. Isso é um alerta honesto: a leitura local ingênua
    simplesmente quebra onde a renovável é dominante. Vale sublinhar o limite do
    exercício — correlação não identifica a curva de oferta nem prova causalidade; o
    valor da figura é apenas documentar o canal que o clearing, adiante, tornará
    estrutural em vez de imposto.
  ]
)

A figura anterior mostrou que o preço reage à carga residual, mas com muito ruído.
A próxima serve para diagnosticar de onde vem esse ruído: queremos saber se o preço
tem uma forma horária estável, que justificaria trabalhar com uma curva média por
hora, ou se ele é tão disperso dentro de cada hora que uma média esconde mais do que
revela. A resposta orienta uma decisão metodológica concreta — se a dispersão for
grande, não podemos representar o ano por um único perfil de preço, e a compressão
por medoides sazonais passa a ser necessária, não conveniente.

#imgfig(
  [Perfil horário do CMO observado com média e dispersão.],
  "outputs/figures/rel_cmo_perfil_horario_erro.png",
  note: [
    O ponto central é a média do CMO em cada hora do dia e a barra cobre mais ou
    menos um desvio-padrão, calculado separadamente por subsistema. O veredito é
    desconfortável para quem gostaria de uma curva de preço simples: a dispersão é
    da ordem da própria média, e é maior justamente nas horas críticas do entardecer,
    onde as decisões de flexibilidade se concentram. Na prática, isso significa que
    duas horas "iguais" no relógio podem ter preços muito diferentes conforme a
    hidrologia e a carga do dia, e é essa a razão de fundo para não colapsar 2025 num
    perfil médio único e para levar a variabilidade sazonal para dentro dos medoides.
    Uma ressalva técnica que evita leitura excessiva: as barras medem variabilidade
    histórica, não são intervalos de confiança e não pressupõem normalidade; a cauda
    de preços altos é assimétrica e a barra simétrica a subestima.
  ]
)

Até aqui tratamos nível e rampa residual como conceitos. A curva do pato é a figura
que os torna visíveis de uma vez só e conecta os dois: queremos ver, no Sudeste/CO,
como a solar rebaixa o nível no meio do dia e, ao mesmo tempo, íngrime a rampa do
entardecer. É a tradução gráfica da distinção da abertura desta seção — capacidade
versus flexibilidade — e é o fato estilizado que motiva boa parte do resto do
relatório: não basta ter energia, é preciso tê-la na hora certa.

#imgfig(
  [Curva do pato no Sudeste/CO: carga, líquida de solar e residual.],
  "outputs/figures/rel_curva_pato.png",
  note: [
    A curva preta é a carga média ao longo do dia; a laranja subtrai a solar,
    revelando a carga líquida; a verde subtrai também a eólica, chegando à carga
    residual. O vale residual médio afunda para cerca de 30,1 GW por volta das 12h —
    é a solar "escavando" o meio do dia. O ponto que realmente importa vem depois:
    quando o sol se põe, a maior subida média de uma hora chega a cerca de 5,17 GW
    entre 16h e 17h, e essa rampa tem de ser coberta em uma hora por hidro e térmica.
    Ou seja, a solar não elimina o problema, ela o *desloca* — troca escassez de
    energia diurna por necessidade de flexibilidade no entardecer, exatamente o
    trade-off que a modelagem precisa capturar. Duas cautelas de leitura são
    importantes aqui: primeiro, este é um perfil *médio*, e a rampa em dias
    individuais ruins é bem mais severa do que 5,17 GW — por isso o clearing trabalha
    com os medoides, não com a média; segundo, a figura é do Sudeste/CO, onde a solar
    domina, e não se transporta para o Nordeste, cujo formato é governado pela eólica.
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

O mesmo fenômeno pode ser derivado, e não apenas medido, o que dá confiança de que
não é um artefato da amostra. Escreva o preço como função da carga residual,
$lambda_t = Lambda(L_t - K^"agg" a_t^s - G_t^w)$, com $Lambda' > 0$. Aqui $Lambda$
(lambda maiúsculo) é o *mapa* que traduz carga residual em preço — a versão
estilizada do clearing —, e $lambda_t$ (lambda minúsculo), em R\$/MWh, é o *valor*
que esse mapa assume na hora $t$; a hipótese $Lambda'>0$ apenas formaliza o fato
empírico de que menos carga residual significa preço mais baixo. O termo
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
    A linha de referência em um separa quem captura acima da média do sistema de quem
    captura abaixo. O resultado é nítido e na direção prevista: a solar fica abaixo
    de um nos quatro subsistemas, entre 0,640 e 0,729, enquanto a eólica fica acima,
    entre 1,044 e 1,309. Em números, a solar do Nordeste recebe cerca de 36% menos do
    que a média do sistema — um desconto grande, que já existe *sem* a expansão que o
    modelo vai estudar, sinal de que o ponto de partida já é desfavorável e tende a
    piorar. Vale notar o contraste locacional: o desconto solar é mais fundo no
    Nordeste (0,640) do que no Sudeste/CO (0,725), o oposto do que a distribuição de
    capacidade sugeriria, e a eólica nordestina é a que mais captura acima da média —
    a complementaridade temporal entre vento e carga trabalha a favor dela na
    amostra. A cautela central é que isto é um fator calculado com CMO, não um retorno
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
chamam exatamente este mesmo operador. Nenhuma delas troca o despacho por baixo dos
panos; elas diferem apenas em *quem escolhe a capacidade* e em *como essa escolha
volta a alimentar o campo agregado*.

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

== Rede e variáveis de decisão

Há quatro nós físicos e um hub contábil $I$. As arestas são

$ I arrow.l.r "SE", quad I arrow.l.r "NE", quad I arrow.l.r "N",
  quad "SE" arrow.l.r "S". $

O hub não possui carga nem geração e satisfaz conservação algébrica dos fluxos. Os
limites reduzidos são 8,0 GW em I-SE, 8,0 GW em I-NE, 4,5 GW em I-N
e 6,0 GW em SE-S. Não existem ligações diretas N-NE ou N-SE no
grafo reduzido.

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
    limite bilateral — não dois fluxos independentes. A leitura estrutural importa
    mais do que parece: Norte, Nordeste e Sudeste/CO só se comunicam *através* do
    hub, e o Sul pendura-se exclusivamente no Sudeste/CO. Isso tem consequência
    direta — não há caminho N–NE ou N–SE direto, então o excedente eólico do
    Nordeste precisa passar pelo hub para escoar, e é aí que os 8,0 GW da interface
    I–NE viram um gargalo que pode "trancar" energia barata na origem. Essa é a
    representação mais frágil do modelo e convém dizê-lo sem rodeios: quatro
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
    energia equivalente agregada de todo o subsistema, e um traço bonito aqui não
    valida a hidrologia de nenhuma usina específica.
  ]
)

== Problema de despacho

O que o clearing minimiza é o *custo total de operar o sistema* ao longo do ano, e
cada parcela da função objetivo é um custo econômico ou uma penalidade que
representa uma preferência operacional. Escrevendo $p_d$ para o peso do medoide (os
dias que ele representa, o que anualiza o custo),

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

Lendo termo a termo: $c_(ell,b)^"th" n_(d,ell,b,t)$ é o custo do combustível
térmico, o preço $c_b^"th"$ (R\$/MWh) do bloco vezes o quanto ele gera $n_b$ (MWh) —
é o termo que ordena o mérito; a pequena parcela $frac(epsilon_n,2) n^2$ apenas
convexifica cada bloco. O termo $pi_h h$ cobra um custo de oportunidade $pi_h$ por
turbinar água agora em vez de guardá-la. As parcelas $pi_s c^s$ e $pi_w c^w$
penalizam levemente o corte solar e eólico, e $pi_u u$ penaliza *pesadamente* o
déficit $u$ — a energia não atendida — com um valor próximo ao custo social de não
atender a carga. O termo $gamma_r (n_t - n_(t-1))^2$ desincentiva variações bruscas
da térmica entre horas, representando o custo de rampa. Por fim, $epsilon_x ||x||_2^2$
é um regularizador numérico sobre o vetor $x$ de todas as variáveis de despacho, que
existe para tornar a solução única e estável.

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
custo marginal. Não há troca posterior de sinal. Em um bloco térmico interior e omitindo apenas
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
    30% de afluência não encarece o sistema em 30%, encarece-o em quase vinte vezes —
    o preço médio salta de 321 para 5.834,95 R\$/MWh e surgem 71,69 TWh de déficit —
    enquanto adicionar 30% mal derruba o preço, de 321 para 241. Essa é a marca de um
    sistema que opera perto do limite de adequação: há muito a perder na seca e pouco
    a ganhar na cheia. Vale uma leitura fina de dois números que parecem contradizer
    a intuição. Primeiro, os 5.834,95 R\$/MWh ainda ficam *abaixo* do custo de
    déficit de 8.327 — não porque falte escassez, mas porque é uma média sobre nós e
    horas, e muitas horas seguem baratas; a média esconde os picos onde o preço bate
    o teto. Segundo, no cenário úmido o vertimento dispara para 130,47 TWh: água em
    excesso não vira valor, vira desperdício, porque o reservatório equivalente e o
    fechamento cíclico não deixam guardá-la para depois. As ressalvas pesam: o choque
    é uniforme e determinístico, não representa uma distribuição de energia natural
    afluente nem a dependência entre bacias, e uma seca real, concentrada em algumas
    regiões e persistente por meses, teria efeito distinto do que este corte parelho
    de 30% sugere.
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
clearing, a receita anual do produtor (em R\$) e o preço que ele captura (em
R\$/MWh) são

$ Rev_i=sum_d p_d sum_t lambda_(d,"SE",t) g_(i,d,t)^s, $

$ P_i^"cap"=frac(Rev_i, sum_d p_d sum_t g_(i,d,t)^s). $

O lucro desconta o custo anualizado da capacidade, $q_("SE",theta_i)$ em R\$/MW/ano
para o tipo $theta_i$ do projeto:

$ Pi_i=Rev_i-q_("SE",theta_i)K_i. $

Um detalhe de honestidade numérica: quando $kappa=0$ o entrante simplesmente não
existe, então sua receita é zero e seu preço capturado é *indefinido* — o
denominador também é zero. Preferimos deixar a lacuna a preenchê-la com um zero
artificial que distorceria a leitura da curva.

Queremos ver, na figura, duas coisas ao mesmo tempo: se a entrada corrói a receita
do incumbente (a canibalização) e se essa corrosão é gradual ou tem saltos. A
resposta ao segundo ponto é a mais reveladora e não é óbvia à primeira vista.

#imgfig(
  [Incumbente e entrante: preço capturado, lucro e curtailment.],
  "outputs/figures/rel_finite_two_player.png",
  note: [
    No eixo horizontal está $kappa=K_2/K_1$; os painéis mostram, em sequência, o
    preço capturado, o lucro anual e o corte alocado a cada projeto. A história tem
    dois regimes muito diferentes, e é isso que importa. Enquanto o sistema ainda tem
    déficit (de $kappa=0$ a $0,50$), acrescentar solar quase não mexe no preço
    capturado do incumbente — cai só de 3.011 para 2.967 R\$/MWh — porque cada MWh
    novo está substituindo energia não atendida, que é caríssima; nesse trecho o
    déficit despenca de 1,615 TWh para 0,00183 TWh. Então, entre $kappa=0,50$ e
    $0,75$, o sistema cruza o *limiar de adequação*: o déficit some e, de repente, o
    preço capturado desaba para 369 R\$/MWh, um oitavo do que era. É um penhasco, não
    uma ladeira. A lição econômica é forte: o valor da solar não é uma propriedade
    fixa da tecnologia, ele depende inteiramente de o sistema estar ou não com aperto
    — os primeiros MW valem ouro, os seguintes valem uma fração. A cautela é igualmente
    importante e limita o que se pode concluir: essa descontinuidade é um efeito de
    *mudança de restrição ativa* (a saída do déficit), não a prova de uma curva de
    preço suave, e como aqui os dois produtores são grandes e movem o preço de
    propósito, isto ainda *não* é o campo médio — é a intuição que ele vai
    disciplinar. O curtailment, aliás, permanece desprezível em toda a faixa, então o
    salto não vem de corte, vem de preço.
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
$a^s$ vezes a produtividade $xi_theta$ — e pela fração de fato entregue
$delta_(d,ell,t)^"delivery" = g^s slash overline(g)^s$, que desconta o curtailment.
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
$kappa_u$ e multiplicado pela capacidade $k$. O segundo termo, $frac(gamma,2)k^2$, é
uma penalidade quadrática crescente que representa custos marginais de expansão que
sobem com o tamanho: $gamma$ controla essa curvatura e é, na prática, o que impede
que o investidor queira capacidade infinita quando o lucro por MW é positivo.
Maximizando, a capacidade ótima interior é

$ k^* = pos(frac((R-q+p)kappa_u,gamma)), $

onde $pos(dot)$ zera valores negativos (não se investe quando o lucro por MW é
negativo), truncada ao intervalo da grade. Quando $k^*$ cai entre dois nós $b$ e
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
    alvo de calibração. A rotina encerra na iteração 146 com
    $r_C=5,25 times 10^(-5)$, $r_m=4,75 times 10^(-4)$ e $r_R=3,65 times 10^(-7)$, os
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
    azul é o planejador social, que escolhe 77,91 GW — praticamente o *dobro*. O
    resultado é grande demais para ser lido ingenuamente, e três leituras críticas
    são obrigatórias. Primeiro, a comparação é honesta porque os três primeiros são
    ancorados no mesmo total, então o salto do planejador não é artefato de escala,
    é o valor que ele atribui à adequação evitada. Segundo, e decisivo: o Sudeste/CO
    social bate *exatamente* em 57,23 GW, que é três vezes a proxy — ou seja, o teto
    regional $3K^"obs"$ está ativo. Isso muda a interpretação: 77,91 GW é o ótimo de
    um problema *limitado*, e o planejador irrestrito quereria ainda mais no
    Sudeste/CO; o número não é um ponto interior bem-comportado, é uma quina imposta
    pela restrição que escolhemos. Terceiro, este resultado *inverte* a conclusão
    preliminar do relatório antigo, o que obriga cautela redobrada — ele diz que,
    nesta parametrização, o benefício social de mais solar (sobretudo por adequação)
    supera seu custo bem além do que o privado investe. Mas "nesta parametrização"
    carrega muito: o resultado é sensível ao custo de déficit, ao custo anualizado,
    ao teto e à ausência de contratos. A leitura correta não é "devemos dobrar a
    capacidade", e sim "há um fosso entre valor privado e social grande o bastante
    para justificar investigar por quê".
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
degrau vale 93 MW. A *ação* combina para onde ir e quanto crescer,

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

Como $m_y$ é uma fração, precisamos de uma escala para traduzi-la em MW. Chamando $M$
o número equivalente de agentes, a capacidade instalada em cada nó é

$ K_(ell,y)(m_y)=M sum_(theta,b) K_b m_y(ell,theta,b),
  quad K_b=93b " MW", $

ou seja, a soma da capacidade $K_b$ de cada bin ponderada pela fração de agentes ali.
O valor de $M$ é fixado uma única vez, para que o estado inicial — concentrado em
$b=4$ — reproduza os 39,069 GW da proxy observada. Isso ancora a escala da população,
mas não impõe nada sobre a trajetória que virá depois.

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
modelo agregado que não corresponde a como projetos reais se distribuem. É uma
muleta modeladora honesta, não um dado.

== Bellman regularizado e política

Um agente que decide hoje precisa pesar o retorno imediato contra o valor de onde a
ação o deixa amanhã. O objeto que faz essa conta é o *valor da ação* $Q_y(x,a)$: a
recompensa imediata mais o valor futuro descontado do estado para onde a transição
leva,

$ Q_y(x,a)=r_y(x,a;m_y,R_y)+beta V_(y+1)(T(x,a)), $

onde $beta$ desconta o futuro e $V_(y+1)$ é o *valor do estado* no estágio seguinte.
Em vez de simplesmente pegar a melhor ação, adotamos uma versão *regularizada*: com
um parâmetro $tau>0$, o valor do estado é uma média suave (o `log-sum-exp`) sobre as
ações viáveis $cal(A)(x)$,

$ V_y(x)=tau log sum_(a in cal(A)(x)) exp(Q_y(x,a)/tau), $

e a política — a probabilidade de escolher cada ação — é o `softmax` correspondente,

$ pi_y(a|x)=frac(exp(Q_y(x,a)/tau),
  sum_(a' in cal(A)(x)) exp(Q_y(x,a')/tau)). $

O papel de $tau$ é o de uma *temperatura*: quanto maior, mais a política se espalha
entre ações de valor parecido; quanto menor, mais ela se concentra na melhor. O
baseline usa $tau=100$ milhões de R\$, um valor alto, então o objeto que calculamos
é explicitamente um *MFG entrópico regularizado*. Isso merece destaque porque muda a
interpretação: a aleatoriedade da política não é ruído acrescentado depois de otimizar
— ela *é* a solução ótima do problema regularizado. Uma consequência que não se deve
varrer para baixo do tapete é que um $tau$ grande suaviza escolhas de propósito, e a
dispersão que veremos nos mapas de política é, em parte, essa regularização, não
necessariamente indiferença econômica real.

Quando $tau=0$, o código substitui o `log-sum-exp` pelo máximo de Bellman e divide
probabilidade apenas entre ações numericamente empatadas:

$ pi_y(a|x)>0 arrow.r Q_y(x,a)=max_(a')Q_y(x,a'). $

Não se afirma que o equilíbrio regularizado com $tau>0$ seja idêntico ao jogo
discreto não regularizado.

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
a comparar com a que a gerou. Como a transição é determinística e cada política soma
um, a equação conserva massa por construção — nenhum agente some ou é criado. Em memória, o notebook verifica numericamente

$ max_y abs(sum_x m_y(x)-1)=2,22 times 10^(-16), $

e, para cada estado e estágio,

$ max_(y,x) abs(sum_a pi_y(a|x)-1)=4,44 times 10^(-16). $

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

mas a certificação é feita no operador sem relaxação:

$ r_m=||hat(m)-m||_infinity, $

$ r_C=max_(y,ell) frac(abs(hat(K)_(y,ell)-K_(y,ell)),K_(y,ell)+1), $

$ r_R=max_(y,ell,theta)
  frac(abs(hat(R)_(y,ell,theta)-R_(y,ell,theta)),
       abs(R_(y,ell,theta))+1). $

As tolerâncias são $2 times 10^(-4)$ para $r_m$ e
$5 times 10^(-3)$ para $r_C$ e $r_R$.

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
    com $r_m=1,91 times 10^(-4)$, $r_C=3,90 times 10^(-4)$ e $r_R=5,85 times 10^(-6)$,
    os três dentro das tolerâncias: existe uma trajetória internamente consistente.
    No plano econômico, a capacidade sobe de 39,069 para 55,376 GW, com a expansão
    fortemente concentrada no Sudeste/CO — o mesmo destino que o planejador e o
    benchmark privado já apontavam, agora com localização e valor futuro escolhidos
    juntos. Duas ressalvas críticas equilibram o resultado. A primeira: 55,376 GW fica
    *abaixo* dos 77,91 GW do planejador social e *acima* dos 39,088 GW do privado
    estático, o que faz sentido — o agente dinâmico enxerga valor futuro que o privado
    míope não vê, mas continua sem internalizar as externalidades que só o planejador
    considera. A segunda, e é a que mais importa: nada disso valida os parâmetros. O
    número terminal é fortemente condicionado por $tau$, pelos custos de entrada e
    realocação e pelo horizonte de seis estágios, e uma escolha diferente desses
    valores moveria o destino; o que a figura certifica é o ponto fixo, não que ele
    seja o do sistema real.
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

A trajetória agregada esconde *como* a expansão acontece — se todos crescem um pouco
ou se alguns saltam para bins altos. A figura seguinte abre a distribuição inteira
para mostrar esse mecanismo, e o motivo de olhá-la é justamente distinguir crescimento
por adesão de crescimento por intensificação.

#imgfig(
  [Trajetória completa da distribuição por localização e capacidade individual.],
  "outputs/figures/rel_dynamic_mfg_distribution_trajectory.png",
  note: [
    Cada painel é um subsistema; o eixo horizontal é o estágio $y=0,dots,6$, o
    vertical é a capacidade individual $K_b$ e a cor é a massa $sum_theta m_y(ell,
    theta,b)$ naquele nó e bin, com escala comum aos quatro. Toda a população começa
    concentrada em 372 MW por agente e, ao longo do horizonte, o Sudeste/CO empurra
    massa para bins de 558, 651, 744 MW e acima — é intensificação, os agentes que já
    estão lá ficam maiores, não apenas mais numerosos. As outras três regiões também
    se espalham, mas com peso agregado bem menor. O ponto que exige cuidado de leitura
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
    localização, integrando tipos e bins. O movimento é claro: o Sudeste/CO ganha
    massa relativa, de 48,83% para 60,08%, enquanto Nordeste e Sul encolhem e o Norte
    cresce de leve. Cruzando com a tabela agregada, aparece a distinção que justifica
    reportar massa e capacidade separadamente: o Nordeste *perde* participação de
    população mas sua capacidade quase não cai (11,789 para 10,978 GW), porque os
    agentes que ficam sobem de bin — menos projetos, cada um maior. A concentração no
    Sudeste/CO é o fio condutor que atravessa todas as camadas do relatório, e aqui
    ela é escolha endógena, não imposição. Um alerta de interpretação para não ler
    demais na figura: isto é realocação de *agentes de investimento*, não fluxo físico
    de energia entre subsistemas — a rede da Seção 4 é que move elétrons, esta figura
    move decisões.
  ]
)

== Política de localização e capacidade

As figuras anteriores mostraram o *resultado* do movimento; as duas próximas mostram
a *regra de decisão* que o produz — a política em si. Queremos ver de onde vem a
migração para o Sudeste/CO: é atração para lá, ou fuga das outras regiões? E onde a
política decide crescer versus encolher?

#imgfig(
  [Política no estágio inicial: realocação e ajuste esperado.],
  "outputs/figures/rel_dynamic_mfg_policy.png",
  note: [
    As colunas separam os tipos `high` e `low`; na linha de cima a cor é a
    probabilidade de mudar de localização ($ell' != ell$) e na de baixo é o ajuste
    esperado de capacidade $EE[Delta K|x]$ em MW, vermelho para expansão e azul para
    redução, com localizações atuais nas linhas e bins atuais nas colunas. A leitura
    econômica é coerente com a migração observada: estados do Sudeste/CO quase não
    querem sair (probabilidade de realocação perto de zero) e ainda expandem em quase
    toda a grade — é um destino que retém e engorda —, enquanto Nordeste e Sul mostram
    propensão a migrar justamente nos bins altos, onde a canibalização local já
    corrói a receita. Ou seja, o movimento é as duas coisas: atração do Sudeste/CO e
    fuga dos bins saturados das outras regiões. Há um artefato de borda que vale
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
    separando tipos. Como narrativa, ela é útil: dá para seguir um agente
    representativo do Sudeste/CO subindo de bin estágio a estágio, tornando concreto o
    que os mapas agregados dizem em massa. Mas o limite é fundamental e precisa ser
    dito sem meias palavras: a ação modal *joga fora* toda a probabilidade das ações
    não modais, e somar essas trajetórias representativas **não** reconstrói a
    distribuição $m_y$ — uma população em que 51% sobe e 49% desce apareceria aqui
    como "todos sobem", o que é falso no agregado. Qualquer conclusão sobre a
    população tem de vir da equação forward e dos mapas de distribuição; esta figura é
    ilustração pedagógica, não evidência quantitativa.
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
