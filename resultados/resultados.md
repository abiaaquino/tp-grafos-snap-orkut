# Resultados consolidados — G27 com-Orkut

Máquina: Intel(R) Xeon(R) Processor @ 2.10GHz (2 núcleos lógicos), Python 3.13.18, Linux-6.18.44-fc-v726-x86_64-with-glibc2.39.

## Parte I — estrutura

- |V| = 9,917, |E| = 173,413; grau mín./máx./médio/mediano = 1/250/34.97/26; densidade = 0.00353
- 1 componente(s); diâmetro = 9, raio = 5, caminho médio = 3.764 (exato, 98,336,972 pares ordenados)
- triângulos = 1,158,975; clusterização média = 0.3703; transitividade = 0.3276

## Parte II — algoritmos

- Euler: 4,930 vértices de grau ímpar → não euleriano nem semi-euleriano
- Tarjan: 138 pontes e 136 pontos de articulação (confere com networkx: True)
- Dijkstra com peso unitário = BFS: True; Bellman-Ford = Dijkstra: True (8 rodadas de 9,916 possíveis)
- AGM: Prim = Kruskal = 1183.4867 (igraph confere: True)
- Floyd-Warshall: validado (numpy = Python = Dijkstra) em n = 150; no grafo completo exigiria 750 MB e ~0.4 h (numpy) ou ~0.33 dias (Python puro)

| Algoritmo | Complexidade | média (ms) | DP (ms) | IC95 (ms) | método |
|---|---|---|---|---|---|
| BFS | O(n+m) | 14.545 | 2.606 | [13.613; 15.478] | Normal (z) |
| DFS | O(n+m) | 13.764 | 1.233 | [13.323; 14.206] | Normal (z) |
| Euler (verificação) | O(n+m) | 14.448 | 0.708 | [14.194; 14.701] | Normal (z) |
| Dijkstra | O((n+m) log n) | 62.274 | 7.316 | [59.656; 64.892] | Normal (z) |
| Bellman-Ford | O(n·m) [efetivo O(k·m)] | 161.381 | 10.044 | [157.787; 164.975] | Normal (z) |
| Tarjan (pontes/articulações) | O(n+m) | 31.971 | 1.577 | [31.407; 32.536] | Normal (z) |
| Prim | O(m log n) | 259.642 | 7.970 | [256.790; 262.495] | Normal (z) |
| Kruskal | O(m log m) | 125.872 | 11.431 | [121.781; 129.962] | Normal (z) |

## Parte III — redes complexas

- Small-world vs ER: C/C_rand = 105.1, L/L_rand = 1.300, σ = 80.9, ω = 0.261 → small-world: True
- Small-world vs configuração (graus preservados): C/C_conf = 34.39, L/L_conf = 1.283, σ = 26.81
- Lei de potência (MLE): subgrafo γ = 6.00 (k ≥ 119, cauda 259), vs exponencial: exponencial (p = 0.00358); Orkut completo γ = 2.96 (k ≥ 168, cauda 289,226), vs exponencial: lei de potência (p = 6.08e-252)
- Robustez (r = 496): S aleatória = 0.9992 [0.9991; 0.9993]; ataque Grau: S = 0.9990, c = 8, L = 4.025; ataque Betweenness: S = 0.9960, c = 36, L = 4.246; ataque PageRank: S = 0.9976, c = 22, L = 4.147
- Descoberta: 87.4% das arestas ficam dentro de alguma comunidade comum; C(k) ~ k^-0.20; Spearman(nº comunidades, betweenness) = 0.236

### Modelos clássicos

| Modelo | m | grau máx. | C | L | γ MLE | cauda (Vuong) |
|---|---|---|---|---|---|---|
| Orkut G27 | 173,413 | 250 | 0.3703 | 3.764 | 6.00 | exponencial |
| Erdős–Rényi | 173,413 | 59 | 0.0035 | 2.896 | 19.77 | exponencial |
| Barabási–Albert (m=17) | 168,436 | 568 | 0.0171 | 2.807 | 2.93 | potência |
| Watts–Strogatz (k=34, p=0,1) | 168,589 | 41 | 0.5307 | 3.597 | 20.00 | exponencial |

## Descoberta principal

**Quem sustenta esta região do Orkut não são os hubs, e sim os vértices-ponte entre comunidades.** (i) A clusterização é 34× maior que a de um grafo aleatório com a MESMA sequência de graus, e 87.4% das arestas ligam membros de uma mesma comunidade: o agrupamento vem das comunidades, não dos graus. (ii) Remover os 496 vértices de maior grau (5%) deixa S = 0.9990, praticamente igual à remoção aleatória (S = 0.9992), porque os hubs ficam no núcleo denso e redundante das comunidades. (iii) Remover os 496 de maior betweenness gera 36 componentes (aleatória: 7.8) e aumenta L em 12% (grau: 7%); vértices em 2+ comunidades têm betweenness média 2.2× maior. (iv) A lei de potência existe no Orkut completo (γ = 2.96, k ≥ 168) mas some no subgrafo (cauda exponencial, γ aparente = 6.0): a amostragem por comunidades preserva clusterização e distâncias curtas, mas corta a cauda de graus — os hubs perdem as arestas que saem da amostra.
