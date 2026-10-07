# forestatrisk — referência das funções usadas no "Get started"

Assinaturas extraídas do código-fonte do pacote (v1.3.2). `import forestatrisk as far`.

## Utilitários

| Função | Assinatura | Devolve |
|---|---|---|
| `far.make_dir` | `make_dir(newdir)` | Cria a pasta se não existir |
| `far.countpix` | `countpix(input_raster, value=1, blk_rows=0)` | `dict` com `npix` e `area` (ha) dos pixels = `value` |
| `far.invlogit` | `invlogit(x)` | Logit inverso |

## Dados

```python
far.sample(nsamp=10000, adapt=True, seed=1234, csize=10, var_dir="data",
           input_forest_raster="forest.tif", output_file="sample.txt",
           blk_rows=0, verbose=True)
```
- `input_forest_raster` é o **nome do ficheiro dentro de `var_dir`** (não caminho completo).
- Devolve `pandas.DataFrame`: uma coluna por `.tif` de `var_dir` + `X`, `Y`, `cell`.
- `adapt=True` → nsamp = 1000 por Mha de floresta, limitado a [10 000; 50 000].

Outras funções de dados (pipelines por país via Google Earth Engine):
`far.data.country_download(...)`, `far.data.country_compute(...)`,
`far.check_fcc`, `far.check_aoi`, `far.get_vector_extent`, `far.extent_shp`.

## Modelo

```python
far.cellneigh(raster=None, region=None, csize=10, rank=1)
# -> (nneigh, adj): nº de vizinhos por célula e lista plana de vizinhos (rank=1 → 8 vizinhos)

far.model_binomial_iCAR(
    suitability_formula, data, n_neighbors, neighbors,
    NA_action="drop", data_pred=None, eval_env=0,
    burnin=1000, mcmc=1000, thin=1,
    beta_start=0,        # -99 → estimativas de GLM logístico clássico
    Vrho_start=1,
    mubeta=0, Vbeta=1000,  # priors normais dos betas
    priorVrho=-1.0,      # -1 → 1/Gamma(shape, rate); >0 → valor fixo
    shape=0.5, rate=0.0005, Vrho_max=10,
    seed=1234, verbose=1, save_rho=0, save_p=0)
```
Atributos do objecto: `betas`, `Vrho`, `deviance`, `rho` (média posterior por
célula), `theta_pred`, `mcmc` (array de amostras), `data`, `suitability_formula`.
Métodos: `predict(new_data=None)`, `plot(output_file="mcmc.pdf", plots_per_page=5, ...)`,
`str(model)` → tabela de posteriores (Mean, Std, CI_low, CI_high).

Fórmula: `"I(1-fcc23) + trial ~ scale(v1) + ... + C(pa) + cell"` — `cell` sempre no fim.

Alternativa não espacial: `far.model_random_forest(...)` (scikit-learn).

## Predição

```python
far.interpolate_rho(rho, input_raster, output_file="output/rho.tif",
                    csize_orig=10, csize_new=1)
# escreve também <dir(output_file)>/rho_orig.tif

far.predict_raster_binomial_iCAR(model, var_dir="data",
    input_cell_raster="output/rho.tif", input_forest_raster="data/forest.tif",
    output_file="output/pred_binomial_iCAR.tif", blk_rows=128, verbose=True)
# raster UInt16, probabilidade reescalada 1–65535, 0 = nodata
```
Outras: `far.predict_raster` (modelos sem iCAR), `far.defrate_per_cat`.

## Projecção

```python
stats = far.deforest(input_raster, hectares, output_file="output/fcc.tif", blk_rows=128)
# stats: counts, hectares, threshold, error, error_perc, ndp (nº pixels desmatados), nfp (nº pixels floresta)
```
Outras: `far.allocate_deforestation`, `far.deforest_diffusion`, `far.emissions`
(emissões de carbono a partir de um raster de biomassa/carbono).

## Validação

| Função | Uso |
|---|---|
| `far.computeAUC(pos_scores, neg_scores, n_sample=100000)` | AUC |
| `far.accuracy_indices(pred, obs)` | OA, FOM, Sen, Spe, TSS, K |
| `far.cross_validation(data, formula, mod_type="icar", ratio=30, nrep=5, seed=1234, icar_args={...})` | Validação cruzada |
| `far.map_validation(pred, obs, blk_rows=128)` | Índices de exactidão mapa previsto vs. observado |
| `far.map_confmat`, `far.map_accuracy` | Matriz de confusão por mapas |
| `far.validation_udef_arp`, `far.validation_npix` | Validação por janelas / nº de pixels |
| `far.r_diffproj`, `far.mat_diffproj` | Diferença entre duas projecções |

## Figuras (`far.plot`)

| Função | Assinatura principal |
|---|---|
| `fcc123` | `(input_fcc_raster, output_file="fcc123.png", maxpixels=500000, borders=None, zoom=None, col=[laranja, vermelho, verde], figsize, dpi=300, **kwargs)` |
| `fcc` | `(input_fcc_raster, output_file="fcc.png", maxpixels=500000, borders=None, zoom=None, col_for, col_defor, figsize, dpi=300, **kwargs)` |
| `prob` | `(input_prob_raster, output_file="prob.png", maxpixels=500000, borders=None, legend=False, figsize, dpi=300, **kwargs)` |
| `rho` | `(input_rho_raster, output_file="rho.png", borders=None, figsize, dpi=300, **kwargs)` |

`**kwargs` são passados ao desenho de `borders` (ex. `linewidth=0.2`).
Todas devolvem uma `matplotlib.figure.Figure`.

## Resultados de referência (Guadalupe, seed=1234)

- Amostragem: 168 blocos, 99 células de 10 km (9 × 11).
- Posteriores iCAR: Intercept −3.84; altitude −0.50; slope −0.016; dist_defor −2.06;
  dist_edge −6.89; dist_road −0.041; dist_town −0.092; dist_river −0.012;
  Vrho 3.12; Deviance ≈ 1.52e+04.
- Desmatamento anual 2010–2020: **498.375 ha/ano** → 14 951 ha projectados 2020–2050.
