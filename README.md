# TP — Analisando Grafos do SNAP · G27 **com-Orkut**

**Ana Beatriz Figueiredo de Aquino** · MATA53-393366-2026.2 — Teoria dos Grafos · Instituto de Computação, UFBA

- **Matrícula:** 224215935 → soma dos dígitos = 2+2+4+2+1+5+9+3+5 = **33**; `f(M) = (33 × (35 + 1)) mod 129 = 1188 mod 129 = ` **27**
- **Grafo G27:** [com-Orkut](https://snap.stanford.edu/data/com-Orkut.html): rede social não direcionada, **3.072.441 vértices e 117.185.083 arestas**, com comunidades ground-truth (grupos criados pelos usuários).
- *Nota:* na planilha de datasets o mesmo arquivo do SNAP pode aparecer em mais de um índice/categoria. Este trabalho usa os arquivos `com-orkut.ungraph.txt.gz` e `com-orkut.top5000.cmty.txt.gz`.

## Estrutura

```
scripts/
  tratamento_orkut.py       auditoria do bruto (streaming) + redução por comunidades -> data/processed/
  graus_grafo_completo.py   distribuição de graus do Orkut COMPLETO (para a lei de potência)
  verificar_dados.py        confere data/processed contra o SNAP e o resumo (código ≠ 0 se divergir)
  download_dados.py         baixa os brutos do SNAP e confere o tamanho em bytes
  run_all.py                roda as Partes I–III sem o notebook
  gerar_relatorio.py        monta relatorio/ (macros, tabelas, figuras) a partir de resultados/; --compilar gera o PDF
src/
  tratamento.py   estrutura Grafo (lista de adjacência + CSR), carga, pesos, subgrafos BFS
  algoritmos.py   BFS, DFS, componentes, Euler/Hierholzer, Dijkstra, Bellman-Ford, Floyd-Warshall
                  (numpy e Python puro), Tarjan (pontes/articulações), Prim, Kruskal (união-busca)
  metricas.py     graus, distribuição, densidade, componentes, triângulos, clusterização
  distancias.py   diâmetro/raio/caminho médio exatos (BFS de todos os vértices, igraph em blocos)
  centralidade.py PageRank (iteração de potência) e betweenness de Brandes (exato/amostrado)
  modelos.py      G(n,m), trocas duplas (graus preservados), Barabási–Albert, Watts–Strogatz, small-world (σ, ω)
  leipotencia.py  regressões log-log/bins/CCDF, MLE discreta (zeta de Hurwitz) + KS, Vuong, Poisson
  robustez.py     remoção aleatória × ataques (métricas S, c, L, E) e curvas
  medicao.py      cronometragem (aquecimento, gc desligado), IC Normal/t-Student, ajustes c·f(n,m)
  figuras.py      PNG 300 dpi + PDF          relatorio.py  tabelas LaTeX pt-BR, valores.tex, JSON
  analise.py      orquestra as etapas (usado pelo notebook e pelo run_all.py)
tests/            pytest: implementação própria × networkx/igraph em grafos pequenos com semente fixa
notebook.ipynb    só orquestra; executado do zero, sem erros
data/raw/         bruto do SNAP (NÃO versionado, ver .gitignore)
data/processed/   dataset tratado (3 MB) + resumo/log do tratamento
figuras/          todas as figuras (PNG + PDF)
resultados/       tabelas (CSV + .tex), json/ por etapa, resultados.json, resultados.md, valores.tex
CHECKLIST.md      requisitos do enunciado e onde estão
relatorio/        main.tex (artigo sigconf), referencias.bib, tabelas/, figuras/, 224215935.pdf
```

## Como executar

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate      Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
pytest                                    # 24 testes, ~2 s

# 1) Dados (só se for refazer o tratamento; data/processed já está versionado)
python scripts/download_dados.py          # baixa os 2 arquivos do SNAP para data/raw/ e confere o tamanho
python scripts/tratamento_orkut.py        # 1,5–2,5 min; pico de memória medido: 564 MB (leitura em blocos)
python scripts/graus_grafo_completo.py    # ~1 min
python scripts/verificar_dados.py

# 2) Análises (Partes I–III): ~10 min numa CPU de 2 núcleos
jupyter nbconvert --to notebook --execute --inplace notebook.ipynb
# ou
python scripts/run_all.py                 # --rapido = poucas repetições, só para testar o fluxo

# 3) Relatório (exige LaTeX com acmart; ou envie relatorio/ ao Overleaf, template WebMedia)
python scripts/gerar_relatorio.py --compilar   # -> relatorio/224215935.pdf
```

Os tempos publicados em `resultados/` foram medidos em Linux x86_64, Intel Xeon 2,1 GHz (2 núcleos lógicos), Python 3.13 (registrado em `resultados/resultados.md`).

## Tratamento dos dados (resumo)

| Passo | Resultado |
|---|---|
| Auditoria do bruto (streaming em blocos de 5 M linhas) | 3.072.441 vértices e 117.185.083 arestas, **iguais ao SNAP**; 0 auto-laços; 0 duplicatas (inclusive u–v × v–u); toda linha tem u < v |
| Seleção | comunidade-semente sorteada com `seed=42` (#933) + expansão gulosa pela maior sobreposição \|C∩S\|/\|C\|; comunidades com > 2.000 membros ignoradas |
| Subgrafo induzido por 51 comunidades | todas as arestas do Orkut entre os vértices escolhidos |
| Maior componente conexa (LCC) | **9.917 vértices, 173.413 arestas** (o induzido já era conexo) |

Justificativa: o grafo completo inviabiliza as medidas O(n·m) e O(n³) do TP. O subgrafo induzido por comunidades preserva **todas** as arestas entre os vértices escolhidos (triângulos, graus locais e comunidades), ao contrário da amostragem aleatória de arestas. **Viés:** a amostra fica centrada em grupos densos (clusterização alta) e os graus caem (em média o vértice mantém 40% do grau original). Por isso a lei de potência também foi testada no Orkut completo. Detalhes: `data/processed/resumo_tratamento.md` e `log_tratamento.txt`.

**Pesos.** O Orkut não tem pesos. Usamos (a) **peso unitário**, em que Dijkstra coincide com a BFS e serve de validação, e (b) um **peso derivado dos dados**, `w(u,v) = 1/(1 + |N(u)∩N(v)|)` (*embeddedness*): laços com muitos amigos em comum ficam curtos e arestas-ponte ficam longas. Assim Dijkstra, Bellman-Ford, Floyd-Warshall e a AGM têm uma interpretação.

## Requisito → notebook → código → saída

| Requisito | Seção do notebook | Código | Saída |
|---|---|---|---|
| \|V\|, \|E\|, graus mín./máx./médio | Parte I | `metricas.resumo_graus` | `resultados/parte1_medidas.*` |
| Distribuição de graus | Parte I | `metricas.distribuicao_graus` | `figuras/p1_distribuicao_graus.*` |
| Densidade | Parte I | `metricas.densidade` | `parte1_medidas.*` |
| Componentes e tamanhos | Parte I | `algoritmos.componentes_conexas` | `parte1_componentes.*` |
| Diâmetro, raio, caminho médio | Parte I | `distancias.todos_para_todos` | `parte1_medidas.*`, `parte1_distancias.*`, `p1_distribuicao_distancias.*` |
| Triângulos, clusterização | Parte I | `metricas.triangulos_por_vertice`, `clusterizacao` | `parte1_medidas.*` |
| Visualização reduzida | Parte I | `analise._visualizacao` | `figuras/p1_visualizacao_top150.*` |
| Tabela calculado / não calculado | Parte I | `analise.parte1` | `resultados/parte1_calculadas.*` |
| BFS, DFS | Parte II | `algoritmos.bfs`, `dfs` | `parte2_tempos_*.tex`, `json/parte2.json` |
| Eulerianidade | Parte II | `algoritmos.verificar_euler`, `circuito_euleriano` | `json/parte2.json` |
| Dijkstra, Bellman-Ford | Parte II | `algoritmos.dijkstra`, `bellman_ford` | idem + demonstração de ciclo negativo |
| Floyd-Warshall | Parte II | `algoritmos.floyd_warshall_numpy`, `_python` | `parte2_floyd_warshall.*`, `p2_floyd_warshall.*` |
| Tarjan (pontes/articulações) | Parte II | `algoritmos.tarjan_pontes_articulacoes` | `json/parte2.json` |
| Prim, Kruskal | Parte II | `algoritmos.prim`, `kruskal` | idem |
| Complexidade + tempos (média, DP, IC) | Parte II | `medicao.*`, `analise.parte2` | `parte2_tempos_grafo_completo.*`, `parte2_tempos_escala.*`, `parte2_ajustes.*`, `p2_escala_tempos.*` |
| Small-world (implementação própria) | Parte III.1 | `modelos.small_world`, `erdos_renyi_gnm`, `troca_dupla_arestas` | `parte3_small_world.*`, `p3_small_world.*` |
| Lei de potência | Parte III.2 | `leipotencia.analisar` | `parte3_lei_potencia.*`, `p3_lei_potencia.*` |
| Robustez 5% aleatória × central | Parte III.3 | `robustez.*`, `centralidade.*` | `parte3_robustez*.*`, `p3_robustez_5pct.*`, `p3_robustez_curvas.*` |
| Descoberta | Parte III.4 | `analise.descoberta` | `parte3_hubs.*`, `p3_descoberta.*`, `resultados.md` |
| Extra: ER, BA, WS | Parte III.5 | `modelos.*` | `parte3_modelos.*`, `p3_modelos_ccdf.*` |
| Números para o relatório | Consolidação | `analise.consolidar` | `resultados/resultados.json`, `resultados.md`, `valores.tex` |

## Resultados principais

| Medida | Valor |
|---|---|
| Vértices / arestas | 9.917 / 173.413 |
| Grau mín. / máx. / médio / mediano | 1 / 250 / 34,97 / 26 |
| Densidade | 0,00353 |
| Componentes | 1 |
| Diâmetro / raio / caminho médio (exatos) | 9 / 5 / 3,764 |
| Triângulos / clusterização média / transitividade | 1.158.975 / 0,370 / 0,328 |
| Euler | 4.930 vértices de grau ímpar → nem euleriano nem semi-euleriano |
| Tarjan | 138 pontes (134 levam a folhas) e 136 pontos de articulação |
| Bellman-Ford | converge em 8 rodadas (de 9.916 possíveis); = Dijkstra |
| AGM (peso embeddedness) | Prim = Kruskal = igraph = 1.183,49 |
| Floyd-Warshall no grafo inteiro | 787 MB de matriz; ~0,5 h (numpy) ou ~7,7 h (Python puro), por extrapolação c·n³ (R² > 0,998) |
| Small-world | C/C_ER ≈ 105, L/L_ER ≈ 1,30, σ ≈ 81; contra o modelo de configuração: C/C_conf ≈ 34, σ ≈ 27 → **small-world** |
| Lei de potência | Orkut completo: γ_MLE ≈ 2,96 (k ≥ 168), favorecida contra a exponencial (Vuong); subgrafo: cauda exponencial (corte da amostragem) |
| Robustez (r = 496) | muito robusta: aleatória desprende ~7 vértices da LCC (S ≈ 0,999); ataque por grau desprende 9 (S = 0,999); **ataque por betweenness desprende 38 (c = 36) e alonga L em 12%** |

**Descoberta mais interessante:** os hubs **não** são os pontos críticos desta região do Orkut. Eles ficam no núcleo redundante onde as comunidades se sobrepõem: 73% dos 496 vértices de maior grau pertencem a 2 ou mais comunidades, contra 26% na rede toda. Por isso, removê-los tem quase o mesmo efeito sobre a conectividade que remover vértices ao acaso. Os vértices de maior betweenness são, em sua maioria, outros (só 147 dos 496 coincidem com o top de grau), e removê-los é o ataque mais danoso, ainda que pequeno em termos absolutos. A clusterização (34× a de um aleatório com os mesmos graus) e o fato de 87% das arestas serem intracomunidade mostram que a estrutura vem das comunidades. Os números finais estão em `resultados/resultados.md`.

## Bibliotecas prontas: o que foi usado e para quê

| Biblioteca | Uso |
|---|---|
| numpy / scipy | vetores, CSR, Floyd-Warshall vetorizado; `scipy.special.zeta` (Hurwitz), `scipy.stats` (IC, regressões, Poisson) |
| pandas | leitura em blocos do bruto e tabelas |
| igraph | **apenas** BFS em C para distâncias todos-para-todos (diâmetro, raio, caminho médio, L dos modelos nulos e robustez) e conferência (betweenness exato, PageRank, AGM, transitividade) |
| networkx | montar o grafo da visualização (layout Kamada-Kawai) e servir de gabarito nos testes / conferência de Tarjan |
| matplotlib | figuras |

Todos os algoritmos das Partes II e III (BFS, DFS, Euler, Dijkstra, Bellman-Ford, Floyd-Warshall, Tarjan, Prim, Kruskal, triângulos, PageRank, Brandes, geradores ER/BA/WS, trocas duplas, σ/ω de small-world, MLE da lei de potência) são **implementações próprias** em `src/`.
