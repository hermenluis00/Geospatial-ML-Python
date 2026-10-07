#!/usr/bin/env python3
"""forestatrisk "Get started" pipeline, parameterised.

Reproduces the official tutorial (Guadeloupe) end to end:
data -> sampling -> binomial iCAR model -> rho interpolation ->
deforestation probability -> future forest cover -> figures.

Usage (from an empty working directory):
    python get_started.py                      # Guadeloupe tutorial data
    python get_started.py --data-dir mydata --no-download --csize 20 \
        --period-years 10 --horizon-years 30 --tag 2050

The data directory must follow the tutorial layout (see SKILL.md).
"""

import argparse
import os
import urllib.request
from contextlib import contextmanager
from shutil import copy2
from zipfile import ZipFile

import matplotlib

matplotlib.use("Agg")

import forestatrisk as far  # noqa: E402

DATA_URL = ("https://github.com/ghislainv/forestatrisk/raw/master/"
            "docsrc/notebooks/data_GLP.zip")

DEFAULT_VARIABLES = ["scale(altitude)", "scale(slope)",
                     "scale(dist_defor)", "scale(dist_edge)",
                     "scale(dist_road)", "scale(dist_town)",
                     "scale(dist_river)"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", default="data")
    p.add_argument("--out-dir", default="output")
    p.add_argument("--no-download", action="store_true",
                   help="Do not download the Guadeloupe data set")
    p.add_argument("--nsamp", type=int, default=10000)
    p.add_argument("--csize", type=float, default=10,
                   help="Spatial cell size (km) for the iCAR random effects")
    p.add_argument("--seed", type=int, default=1234)
    p.add_argument("--variables", nargs="+", default=DEFAULT_VARIABLES)
    p.add_argument("--burnin", type=int, default=1000)
    p.add_argument("--mcmc", type=int, default=1000)
    p.add_argument("--thin", type=int, default=1)
    p.add_argument("--period-years", type=float, default=10.0,
                   help="Length of the t2->t3 reference period")
    p.add_argument("--horizon-years", type=float, default=30.0,
                   help="Projection horizon after t3")
    p.add_argument("--tag", default="2050", help="Suffix for projection outputs")
    p.add_argument("--borders", default=None,
                   help="Vector file for map borders (default <data>/aoi_proj.shp or ctry_PROJ.shp)")
    p.add_argument("--blk-rows", type=int, default=10,
                   help="Rows per block for prediction (lower = less memory)")
    p.add_argument("--dpi", type=int, default=800)
    return p.parse_args()


@contextmanager
def forecast_variables(data_dir, names=("dist_edge", "dist_defor")):
    """Temporarily swap dynamic variables for their t3 (forecast) version."""
    swapped = []
    try:
        for n in names:
            cur = os.path.join(data_dir, n + ".tif")
            new = os.path.join(data_dir, "forecast", n + "_forecast.tif")
            if not os.path.exists(new):
                print(f"[warn] {new} not found: using {cur} as is")
                continue
            os.rename(cur, cur + ".bak")
            copy2(new, cur)
            swapped.append(cur)
        yield
    finally:
        for cur in swapped:
            os.remove(cur)
            os.rename(cur + ".bak", cur)


def main():
    a = parse_args()
    d, o = a.data_dir, a.out_dir
    far.make_dir(o)

    # 1. Data
    if not a.no_download:
        if not os.path.exists("data_GLP.zip"):
            urllib.request.urlretrieve(DATA_URL, "data_GLP.zip")
        with ZipFile("data_GLP.zip", "r") as z:
            z.extractall(d)

    dataset = far.sample(nsamp=a.nsamp, adapt=True, seed=a.seed, csize=a.csize,
                         var_dir=d, input_forest_raster="fcc23.tif",
                         output_file=os.path.join(o, "sample.txt"), blk_rows=0)
    dataset = dataset.dropna(axis=0)
    dataset["trial"] = 1
    print(dataset.head(5))

    # 2. Model
    nneigh, adj = far.cellneigh(raster=os.path.join(d, "fcc23.tif"),
                                csize=a.csize, rank=1)
    formula = "I(1-fcc23) + trial ~ " + " + ".join(a.variables) + " + cell"
    mod = far.model_binomial_iCAR(
        suitability_formula=formula, data=dataset,
        n_neighbors=nneigh, neighbors=adj,
        priorVrho=-1, burnin=a.burnin, mcmc=a.mcmc, thin=a.thin,
        beta_start=-99, seed=a.seed)
    print(mod)
    with open(os.path.join(o, "summary_icar.txt"), "w") as f:
        f.write(str(mod))
    mod.plot(output_file=os.path.join(o, "mcmc.pdf"))

    # 3. Predict
    rho_tif = os.path.join(o, "rho.tif")
    far.interpolate_rho(rho=mod.rho, input_raster=os.path.join(d, "fcc23.tif"),
                        output_file=rho_tif, csize_orig=a.csize, csize_new=1)
    prob_tif = os.path.join(o, "prob.tif")
    with forecast_variables(d):
        far.predict_raster_binomial_iCAR(
            mod, var_dir=d, input_cell_raster=rho_tif,
            input_forest_raster=os.path.join(d, "forest", "forest_t3.tif"),
            output_file=prob_tif, blk_rows=a.blk_rows)

    # 4. Project
    fc = [far.countpix(input_raster=os.path.join(d, "forest", f"forest_{t}.tif"),
                       value=1)["area"] for t in ("t2", "t3")]
    with open(os.path.join(o, "forest_cover.txt"), "w") as f:
        f.write("\n".join(str(v) for v in fc) + "\n")
    annual_defor = (fc[0] - fc[1]) / a.period_years
    print(f"Mean annual deforested area t2-t3: {annual_defor} ha/yr")
    fcc_tif = os.path.join(o, f"fcc_{a.tag}.tif")
    stats = far.deforest(input_raster=prob_tif,
                         hectares=annual_defor * a.horizon_years,
                         output_file=fcc_tif, blk_rows=128)
    print({k: v for k, v in stats.items() if k != "counts"})

    # 5. Figures
    borders = a.borders or next(
        (os.path.join(d, f) for f in ("aoi_proj.shp", "ctry_PROJ.shp")
         if os.path.exists(os.path.join(d, f))), None)
    fcc123 = os.path.join(d, "forest", "fcc123.tif")
    if os.path.exists(fcc123):
        far.plot.fcc123(input_fcc_raster=fcc123, maxpixels=1e8,
                        output_file=os.path.join(o, "fcc123.png"),
                        borders=borders, linewidth=0.2,
                        figsize=(5, 4), dpi=a.dpi)
    for name in ("rho_orig", "rho"):
        far.plot.rho(os.path.join(o, name + ".tif"), borders=borders,
                     linewidth=0.5, output_file=os.path.join(o, name + ".png"),
                     figsize=(9, 5), dpi=80)
    far.plot.prob(prob_tif, maxpixels=1e8, borders=borders, linewidth=0.2,
                  legend=True, output_file=os.path.join(o, "prob.png"),
                  figsize=(5, 4), dpi=a.dpi)
    far.plot.fcc(fcc_tif, maxpixels=1e8, borders=borders, linewidth=0.2,
                 output_file=os.path.join(o, f"fcc_{a.tag}.png"),
                 figsize=(5, 4), dpi=a.dpi)


if __name__ == "__main__":
    main()
