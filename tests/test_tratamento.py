import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from tratamento_orkut import chave_canonica, decodificar_chave, selecionar_comunidades


def test_chave_canonica_ignora_sentido_e_detecta_duplicata():
    u = np.array([1, 2, 5, 3]); v = np.array([2, 1, 5, 9])
    k = chave_canonica(u, v)
    assert k[0] == k[1]                      # (1,2) == (2,1)
    assert len(np.unique(k)) == 3
    lo, hi = decodificar_chave(k)
    assert list(lo) == [1, 1, 5, 3] and list(hi) == [2, 2, 5, 9]


def test_selecao_reprodutivel_e_por_sobreposicao():
    cm = [[0, 1, 2, 3], [3, 4, 5], [10, 11], [4, 5, 6, 7], list(range(100, 200))]
    a = selecionar_comunidades(cm, seed=1, alvo_superconjunto=7, tam_max_cmty=50)
    b = selecionar_comunidades(cm, seed=1, alvo_superconjunto=7, tam_max_cmty=50)
    assert a[1] == b[1]                      # mesma semente -> mesma ordem
    assert 4 not in a[1]                     # comunidade grande demais ignorada
    assert len(a[2]) >= 7
