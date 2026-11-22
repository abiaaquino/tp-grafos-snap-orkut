# Checklist do TP — G27 com-Orkut

[R] = precisa aparecer no relatório. Arquivos de saída entre parênteses.

## Dados e tratamento
- [x] [R] Cálculo do grafo pela matrícula (224215935 → 27) — README
- [x] [R] Auditoria do bruto: |V|, |E|, laços, duplicatas, conferência com o SNAP (`data/processed/resumo_tratamento.*`)
- [x] [R] Redução justificada: subgrafo induzido por 51 comunidades ground-truth + LCC (`data/processed/tabela_tratamentos.tex`)
- [x] Bruto fora do Git (`.gitignore`), dados tratados < 50 MB versionados

## Parte I — estrutura
- [x] [R] |V|, |E|, grau mín./máx./médio/mediano, densidade (`resultados/parte1_medidas.*`)
- [x] [R] Distribuição de graus (`figuras/p1_distribuicao_graus.*`)
- [x] [R] Componentes conexas e tamanhos (`resultados/parte1_componentes.*`)
- [x] [R] Diâmetro, raio, caminho médio exatos (`resultados/parte1_distancias.*`, `figuras/p1_distribuicao_distancias.*`)
- [x] [R] Triângulos, clusterização média, transitividade
- [x] [R] Visualização reduzida (`figuras/p1_visualizacao_top150.*`)
- [x] [R] Tabela calculado / não calculado com justificativa (`resultados/parte1_calculadas.*`)

## Parte II — algoritmos
- [x] BFS, DFS, Euler, Dijkstra, Bellman-Ford, Floyd-Warshall, Tarjan (pontes/articulações), Prim, Kruskal (`src/algoritmos.py`)
- [x] Validações cruzadas (Dijkstra × BF × FW; Prim × Kruskal × igraph; Tarjan × networkx)
- [x] [R] Justificativa de pesos (unitário + embeddedness) e de Floyd-Warshall só em subgrafos
- [x] [R] Complexidade teórica + tempos com média, DP, IC (Normal n ≥ 30; t-Student n < 30) (`resultados/parte2_*.tex`)
- [x] [R] Escala em subgrafos e ajuste c·f(n,m) com R² (`figuras/p2_escala_tempos.*`, `figuras/p2_floyd_warshall.*`)

## Parte III — redes complexas
- [x] [R] Small-world com cálculo próprio (σ, ω) vs ER e vs configuração (`resultados/parte3_small_world.*`)
- [x] [R] Lei de potência: regressões, CCDF, MLE + KS, Vuong vs exponencial, Poisson; subgrafo e Orkut completo (`resultados/parte3_lei_potencia.*`)
- [x] [R] Robustez 5% aleatória (30 rep.) × 5% mais centrais (grau, betweenness, PageRank) (`resultados/parte3_robustez*.tex`)
- [x] [R] Descoberta mais interessante (`resultados/parte3_hubs.*`, `figuras/p3_descoberta.*`)
- [x] [R] Extra: ER, BA, WS (`resultados/parte3_modelos.*`, `figuras/p3_modelos_ccdf.*`)

## Entrega
- [x] README com requisito → seção → código → saída
- [x] Notebook executado do zero sem erros
- [x] `resultados/valores.tex` para o relatório
- [x] Relatório `relatorio/224215935.pdf` (falta preencher instituição/e-mail e revisar as interpretações)
- [ ] Commits por etapa e push para o GitHub
