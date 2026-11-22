# Resumo do tratamento — G27 com-Orkut

Gerado por `scripts/tratamento_orkut.py` em 2026-11-22 21:30 (85 s). Log completo em `log_tratamento.txt`.

## Auditoria do grafo original

| Medida | Valor |
|---|---|
| linhas_de_aresta | 117,185,083 |
| vertices_distintos | 3,072,441 |
| id_minimo | 1 |
| id_maximo | 3,072,626 |
| auto_lacos | 0 |
| linhas_com_u_maior_que_v | 0 |
| arestas_nao_direcionadas_unicas_incl_lacos | 117,185,083 |
| linhas_duplicadas_(u,v)_ou_(v,u) | 0 |
| arquivo_ordenado_por_(u,v) | False |
| grau_maximo | 33,313 |
| grau_medio | 76.28142118921079 |
| grau_mediano | 45.0 |
| segundos_leitura | 74.3 |
| linhas_de_comentario | 4 |

Conferência com a página do SNAP: |V| **OK**, |E| **OK**.

## Tratamentos

| Tratamento | Descrição | Resultado |
|---|---|---|
| Leitura do bruto | com-orkut.ungraph.txt.gz, leitura em streaming (blocos de 5,000,000 linhas) | 4 linhas de comentário; 117,185,083 linhas de aresta |
| Auditoria de vértices | vértices distintos que aparecem em alguma aresta | 3,072,441 (SNAP: 3,072,441; ids 1–3,072,626) |
| Auditoria de arestas | linhas de aresta vs. valor oficial | 117,185,083 (SNAP: 117,185,083) |
| Auto-laços | linhas com u = v | 0 |
| Arestas duplicadas | pares repetidos, considerando (u,v) ≡ (v,u) | 0 |
| Orientação | linhas com u > v (arquivo grava cada aresta uma vez?) | 0 |
| Seleção de comunidades | semente aleatória fixa (seed=42) + expansão gulosa por sobreposição |C∩S|/|C|; comunidades com > 2000 nós ignoradas | semente #933; 51 comunidades selecionadas |
| Subgrafo induzido | todas as arestas do grafo original com as duas pontas no conjunto | |V|=9,917; |E|=173,413; 1 componentes |
| Simplificação | remoção de auto-laços e duplicatas; grafo não direcionado, sem pesos | 0 laços; arestas únicas por par |
| Maior componente conexa | mantida só a LCC (distâncias, diâmetro, small-world) | |V|=9,917; |E|=173,413 |

## Grafo final (data/processed)

| Medida | Valor |
|---|---|
| vertices | 9917 |
| arestas | 173413 |
| grau_medio | 34.973 |
| grau_max | 250 |
| grau_min | 1 |
| densidade | 0.0035269135600392494 |
| clusterizacao_media | 0.3703 |
| componentes | 1 |
| auto_lacos | 0 |

## Justificativa metodológica

- O grafo completo (117 M arestas) não cabe confortavelmente em memória com NetworkX e inviabiliza algoritmos O(n²)/O(n³) do TP; por isso foi reduzido.
- Usou-se um **subgrafo induzido por comunidades ground-truth** do próprio SNAP (top5000), e não amostragem aleatória de arestas: o subgrafo induzido preserva **todas** as arestas entre os vértices escolhidos, então graus locais, triângulos e estrutura comunitária são mantidos.
- As comunidades foram encadeadas por **sobreposição** (|C∩S|/|C|) a partir de uma semente aleatória **fixa**, gerando uma região coesa e reprodutível da rede; a LCC garante conectividade para distâncias, diâmetro e small-world.
- **Viés conhecido:** a amostra é centrada em comunidades densas, então clusterização e densidade tendem a ser maiores que no Orkut completo, e o grau no subgrafo é menor que o original (coluna `grau_original` em `orkut_g27_lcc_nos.csv` permite comparar).
- O Orkut é não direcionado e sem pesos: não há versão dirigida a manter; algoritmos ponderados (Dijkstra, AGM) precisarão de pesos declarados explicitamente na análise.

## Arquivos

- `orkut_g27_lcc_arestas.csv` — arestas (src,dst), ids originais do Orkut
- `orkut_g27_lcc.edgelist.gz` — mesma lista em formato edgelist (NetworkX `read_edgelist`)
- `orkut_g27_lcc_nos.csv` — nós, grau original, grau no subgrafo e comunidades
- `orkut_g27_comunidades_selecionadas.csv`, `orkut_g27_prefixos_avaliados.csv` — rastreabilidade da seleção
- `tabela_tratamentos.csv` / `.tex`, `resumo_tratamento.json`, `log_tratamento.txt`
