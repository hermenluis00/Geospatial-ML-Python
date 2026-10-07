---
name: forestatrisk
description: Modelação e projecção espacial do risco de desmatamento com o pacote Python `forestatrisk` (Ghislain Vieilledent, CIRAD), seguindo o tutorial oficial "Get started". Usar sempre que o utilizador pedir para modelar, mapear ou prever desmatamento/desflorestação, calcular a probabilidade espacial de desmatamento, ajustar um modelo logístico binomial com efeitos aleatórios espaciais iCAR, amostrar pixels floresta/desmatado, interpolar efeitos espaciais (rho), projectar a cobertura florestal futura (ex. 2050), calcular a taxa anual de desmatamento, ou produzir mapas fcc123/prob/rho/fcc. Trigger também para "forestatrisk", "far.sample", "model_binomial_iCAR", "deforestation risk map", "risk map REDD+", "jurisdictional baseline / FREL / nível de referência florestal", "mapa de risco de desmatamento", "cenário de desmatamento futuro", ou qualquer pedido de Python + GDAL + desmatamento tropical (Moçambique, Madagáscar, Amazónia, Guadalupe, etc.).
---

# forestatrisk — Modelar e prever o risco de desmatamento

Skill baseada no tutorial oficial **Get started** do `forestatrisk`
(<https://ecology.ghislainv.fr/forestatrisk/notebooks/get_started.html>,
versão do pacote 1.3.x). O caso de estudo é o arquipélago de **Guadalupe**,
com dados de mudança de cobertura florestal 2000–2010–2020.

O fluxo tem 5 etapas: **Dados → Modelo → Predição → Projecção → Figuras**.

```
fcc23.tif + variáveis .tif ──far.sample()──► sample.txt (DataFrame)
                                   │
            far.cellneigh() ──► vizinhança iCAR (nneigh, adj)
                                   │
          far.model_binomial_iCAR() ──► betas, Vrho, rho (por célula 10 km)
                                   │
far.interpolate_rho() ──► rho.tif (1 km)    far.predict_raster_binomial_iCAR() ──► prob.tif
                                   │
far.countpix() ──► taxa anual   far.deforest() ──► fcc_2050.tif
                                   │
                  far.plot.{fcc123, rho, prob, fcc} ──► PNGs
```

## Quando aplicar

- Produzir um **mapa de probabilidade espacial de desmatamento** a partir de
  um raster de mudança florestal e variáveis explicativas.
- Construir um **cenário business-as-usual** (taxa histórica × horizonte)
  e localizar espacialmente o desmatamento futuro (ex. REDD+, FREL,
  alocação jurisdicional).
- Adaptar o tutorial de Guadalupe a outra área de interesse (AOI).

Para o pipeline completo e executável, usar `scripts/get_started.py`.
Para assinaturas detalhadas das funções e objectos devolvidos, ver
`references/api.md`.

## Instalação

Recomendado: ambiente `conda` (o GDAL é a dependência difícil).

```bash
conda create --name conda-far -c conda-forge python gdal numpy matplotlib pandas patsy pip statsmodels earthengine-api --yes
conda activate conda-far
pip install pywdpa scikit-learn   # não disponíveis no conda
pip install forestatrisk          # PyPI
# pip install https://github.com/ghislainv/forestatrisk/archive/master.zip  # versão dev
forestatrisk                      # testa a instalação (imprime o banner)
```

Com `venv`/pip: `python3 -m pip install gdal==$(gdal-config --version)` antes
de `pip install forestatrisk`. Se os bindings GDAL falharem:
`python3 -m pip install --no-cache-dir --force-reinstall gdal==$(gdal-config --version)`.

## Estrutura de dados esperada

```
data/
├── fcc23.tif               # 1 = floresta mantida, 0 = desmatado (t2→t3); nodata fora da floresta
├── altitude.tif  slope.tif
├── dist_defor.tif  dist_edge.tif  dist_road.tif  dist_town.tif  dist_river.tif
├── pa.tif                  # áreas protegidas (0/1)
├── forest/
│   ├── forest_t2.tif  forest_t3.tif   # 1 = floresta
│   └── fcc123.tif                     # mudança 2000–2010–2020 para figura
├── forecast/
│   ├── dist_edge_forecast.tif         # distâncias recalculadas com a floresta em t3
│   └── dist_defor_forecast.tif
└── aoi_proj.shp            # limites (projectados) para os mapas
                            # (no zip actual chama-se ctry_PROJ.shp)
```

Regras críticas:
- **Todos os rasters na mesma projecção métrica, extensão e resolução**
  (o `far.sample` constrói um VRT com todas as `.tif` de `var_dir` como bandas).
- Cada `.tif` em `var_dir` vira uma coluna com o nome do ficheiro — não deixar
  rasters estranhos nessa pasta.
- Os dados de Guadalupe descarregam-se de
  `https://github.com/ghislainv/forestatrisk/raw/master/docsrc/notebooks/data_GLP.zip`.

## Passo a passo

### 1. Dados — amostragem

```python
import os
from shutil import copy2
import forestatrisk as far

far.make_dir("output")

dataset = far.sample(nsamp=10000, adapt=True, seed=1234, csize=10,
                     var_dir="data", input_forest_raster="fcc23.tif",
                     output_file="output/sample.txt", blk_rows=0)
dataset = dataset.dropna(axis=0)   # obrigatório: scale() e o iCAR não aceitam NA
dataset["trial"] = 1               # nº de ensaios binomiais = 1
```

- Amostra **2 × nsamp** pixels (desmatados vs. floresta), estratificados por blocos.
- `adapt=True`: ajusta `nsamp` à área florestal (1000 por Mha, mín. 10 000, máx. 50 000).
- `csize=10`: células espaciais de 10 km para o efeito aleatório iCAR; a coluna
  `cell` identifica a célula de cada ponto. Colunas `X`, `Y` são as coordenadas.

### 2. Modelo — regressão logística binomial com iCAR

```python
nneigh, adj = far.cellneigh(raster="data/fcc23.tif", csize=10, rank=1)

variables = ["scale(altitude)", "scale(slope)",
             "scale(dist_defor)", "scale(dist_edge)", "scale(dist_road)",
             "scale(dist_town)", "scale(dist_river)"]
formula = "I(1-fcc23) + trial ~ " + " + ".join(variables) + " + cell"

mod_binomial_iCAR = far.model_binomial_iCAR(
    suitability_formula=formula, data=dataset,
    n_neighbors=nneigh, neighbors=adj,
    priorVrho=-1,                       # -1 = prior 1/Gamma para Vrho
    burnin=1000, mcmc=1000, thin=1,
    beta_start=-99)                     # -99 = usar estimativas de um GLM simples
print(mod_binomial_iCAR)
```

Pontos-chave:
- Resposta `I(1-fcc23)`: 1 = desmatado, 0 = floresta (inverte o fcc23).
- O termo **`+ cell` tem de ser o último** da fórmula: é interpretado como o
  identificador da célula para o efeito aleatório espacial, não como covariável.
- `csize` em `cellneigh` **tem de ser igual** ao usado em `far.sample`.
- `scale()` (patsy) padroniza as variáveis; os mesmos parâmetros de escala são
  reaplicados automaticamente na predição.
- Variáveis categóricas (ex. `pa`) entram como `C(pa)`.
- Atributos úteis: `.betas`, `.Vrho`, `.rho` (um valor por célula), `.deviance`,
  `.theta_pred` (probabilidade prevista nas observações), `.mcmc` (cadeias),
  `.plot(output_file="output/mcmc.pdf")` (traços MCMC).
- Guardar o resumo: `open("output/summary_icar.txt","w").write(str(mod_binomial_iCAR))`.

Os valores exactos variam ligeiramente entre versões/execuções (MCMC); em
v1.3.2 obteve-se ex. dist_edge −6.99, dist_defor −2.17, Vrho 3.14.

Interpretação típica (Guadalupe): `dist_edge` (−6.89) e `dist_defor` (−2.06)
dominam — o desmatamento concentra-se junto à orla florestal e a desmatamento
recente; `altitude` negativo; `slope`, `dist_road`, `dist_river` com IC a
cruzar 0 (não significativas). Variáveis cujo intervalo credível inclui 0
podem ser removidas e o modelo reajustado.

### 3. Predição

```python
# 3.1 Efeitos aleatórios espaciais: 10 km → 1 km por interpolação cubic spline
far.interpolate_rho(rho=mod_binomial_iCAR.rho, input_raster="data/fcc23.tif",
                    output_file="output/rho.tif", csize_orig=10, csize_new=1)
# (escreve também output/rho_orig.tif com os valores originais por célula)

# 3.2 Actualizar variáveis dinâmicas para a data t3, prever, e repor os originais
os.rename("data/dist_edge.tif", "data/dist_edge.tif.bak")
os.rename("data/dist_defor.tif", "data/dist_defor.tif.bak")
copy2("data/forecast/dist_edge_forecast.tif", "data/dist_edge.tif")
copy2("data/forecast/dist_defor_forecast.tif", "data/dist_defor.tif")

far.predict_raster_binomial_iCAR(
    mod_binomial_iCAR, var_dir="data",
    input_cell_raster="output/rho.tif",
    input_forest_raster="data/forest/forest_t3.tif",
    output_file="output/prob.tif",
    blk_rows=10)                       # blocos pequenos = menos memória

os.remove("data/dist_edge.tif"); os.remove("data/dist_defor.tif")
os.rename("data/dist_edge.tif.bak", "data/dist_edge.tif")
os.rename("data/dist_defor.tif.bak", "data/dist_defor.tif")
```

- A predição usa a floresta em **t3 (2020)**, logo as distâncias dinâmicas
  (`dist_edge`, `dist_defor`) devem refletir t3 — daí a troca de ficheiros.
  Usar `try/finally` para garantir a reposição mesmo em caso de erro.
- `prob.tif` é UInt16 reescalado **1–65535** (0 = nodata/fora da floresta).

### 4. Projecção da cobertura florestal futura

```python
fc = [far.countpix(input_raster=f"data/forest/forest_{d}.tif", value=1)["area"]
      for d in ("t2", "t3")]            # área em ha
annual_defor = (fc[0] - fc[1]) / 10.0   # período 2010–2020 → 498.375 ha/ano em Guadalupe

stats = far.deforest(input_raster="output/prob.tif",
                     hectares=annual_defor * 30,   # 2020–2050
                     output_file="output/fcc_2050.tif", blk_rows=128)
```

- `far.deforest` ordena os pixels pela probabilidade e desmata os de maior
  risco até atingir os hectares pedidos (procura de limiar que minimiza o erro).
- `stats` → `counts`, `hectares`, `threshold`, `error`, `error_perc`, `ndp`, `nfp`.
- Hipótese implícita: taxa constante (cenário BAU). Para outros cenários,
  alterar `hectares` (ex. taxa ajustada, intervalos 2030/2040/2050 em ciclo).

### 5. Figuras

```python
far.plot.fcc123(input_fcc_raster="data/forest/fcc123.tif", maxpixels=1e8,
                output_file="output/fcc123.png", borders="data/aoi_proj.shp",
                linewidth=0.2, figsize=(5, 4), dpi=800)
far.plot.rho("output/rho_orig.tif", borders="data/aoi_proj.shp", linewidth=0.5,
             output_file="output/rho_orig.png", figsize=(9, 5), dpi=80)
far.plot.rho("output/rho.tif", borders="data/aoi_proj.shp", linewidth=0.5,
             output_file="output/rho.png", figsize=(9, 5), dpi=80)
far.plot.prob("output/prob.tif", maxpixels=1e8, borders="data/aoi_proj.shp",
              linewidth=0.2, legend=True, output_file="output/prob.png",
              figsize=(5, 4), dpi=800)
far.plot.fcc("output/fcc_2050.tif", maxpixels=1e8, borders="data/aoi_proj.shp",
             linewidth=0.2, output_file="output/fcc_2050.png",
             figsize=(5, 4), dpi=800)
```

Em ambientes sem display: `import matplotlib; matplotlib.use("Agg")` antes de
importar `forestatrisk`.

## Adaptar a outra área (ex. Moçambique)

1. Preparar `fcc23.tif`, `forest_t2/t3.tif`, `fcc123.tif` (ex. TMF/JRC, GFC
   Hansen) e as variáveis, todas reprojectadas para a mesma CRS métrica
   (ex. UTM 36S, EPSG:32736) e grelha — o próprio `forestatrisk` tem
   `far.data.country_download` / `far.data.country_compute` (via GEE) para isto.
2. Ajustar `csize` à dimensão da AOI (10 km é um bom ponto de partida;
   AOIs muito grandes → 20–30 km para limitar o nº de células).
3. Rever a lista `variables` à disponibilidade local; adicionar `C(pa)` se
   houver áreas protegidas.
4. Para áreas grandes: `blk_rows` pequeno na predição, `adapt=True` na amostragem.
5. Ajustar `T` (anos do período de referência) e o horizonte da projecção.
6. Validar: `far.computeAUC`, `far.accuracy_indices`, `far.cross_validation`,
   e validação do mapa com `far.map_validation` (ver `references/api.md`).

## Erros comuns

| Sintoma | Causa / solução |
|---|---|
| Erro em `scale()` ou no iCAR com NaN | Faltou `dataset.dropna(axis=0)` |
| `ModuleNotFoundError: osgeo` | GDAL Python não instalado ou versão ≠ `gdal-config --version` |
| `_ARRAY_API not found` / `numpy.core.multiarray failed to import` | Bindings GDAL do sistema compilados com NumPy 1.x: reinstalar no venv com `pip install --no-build-isolation "gdal[numpy]==$(gdal-config --version)"` (NumPy 2 já instalado); não fazer downgrade de numpy |
| `RuntimeError: data/aoi_proj.shp: No such file` | O zip actual de Guadalupe traz `data/ctry_PROJ.shp`; usar esse em `borders=` |
| Coluna inesperada no DataFrame | Raster extra em `var_dir`; mover para subpasta |
| Índices de `cell` incoerentes | `csize` diferente entre `sample` e `cellneigh` |
| `MemoryError` na predição | Reduzir `blk_rows` (ex. 10) |
| Mapas de probabilidade irrealistas | Não actualizou `dist_edge`/`dist_defor` para t3 |
| `data/dist_*.tif.bak` deixados para trás | Execução interrompida — repor manualmente os `.bak` |

## Referência

Vieilledent G. (2021). *forestatrisk: a Python package for modelling and
forecasting deforestation in the tropics.* Journal of Open Source Software,
6(59), 2975. <https://doi.org/10.21105/joss.02975>
