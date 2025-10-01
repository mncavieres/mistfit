# Nested Sampling MIST fitter

Nested-sampling stellar-parameter fits using **dynesty** + **minimint** (MIST isochrones with extinction).

This package provides a single public function, `fit_stars_with_minimint`, which takes an `astropy.table.Table` containing photometry (and optionally spectroscopy, parallax, and reddening priors), and returns the same table with posterior summaries appended. A simple CLI is also included.

---
## Citation

If this tool is used in a publication, please cite **dynesty**, **minimint**, and Cavieres (in prep).


## Quick start

### Python example

```python
from astropy.table import Table
from mistfit import fit_stars_with_minimint

if __name__ == "__main__":
    tab = Table.read("xshooter_with_phot_gaia_minimint.fits")
    out = fit_stars_with_minimint(
        tab,
        output_path="xshooter_nested",
        nlive=5000,
        processes=8,
        debug=True,
    )
```

### CLI

```bash
python minimint_fit.py /path/to/input_table.fits \
  --outdir /path/to/out \
  --nlive 5000 \
  --dlogz 0.01 \
  --procs 8 \
  --debug
```

The outputs are written to `--outdir`:

* Updated table with appended posterior columns: `fit_results.fits` (or `.ecsv` if FITS types conflict)
* If `--debug`, per-star folders with `runplot.png`, `traceplot.png`, `corner.png`, optional `isochrone_check.png`, and `summary.json`.

---

## Input table format

**Minimum requirement**: `≥ 3` photometric bands recognized by *minimint* **and** their per-band uncertainties, using the naming conventions below. Each row corresponds to one star. A `source_id` column is recommended (used to name per-star debug folders), but not required.

### Column names and uncertainty suffixes

* **Photometric band columns** use the *minimint* band names (see list below). Examples: `Gaia_G_EDR3`, `DECam_g`, `2MASS_Ks`.
* The **uncertainty column** for a quantity `<NAME>` is discovered by trying these suffixes (first match wins):

  * `<NAME>_ERR`, `<NAME>_ERRMAG`, `<NAME>_SIG`, `<NAME>_E`

  Examples: `DECam_g_ERR`, `Gaia_BP_EDR3_ERR`, `Teff_ERR`, `PARALLAX_E`.

### Supported photometric bands (subset)

The code auto-detects any of the following present in your table (case-sensitive):
`Bessell_I, PS_z, 2MASS_Ks, SkyMapper_u, Gaia_RP_MAW, Gaia_RP_DR2Rev, Gaia_BP_EDR3, PS_w, PS_r, SDSS_g, SkyMapper_v, Bessell_V, PS_y, SDSS_r, Bessell_U, Gaia_BP_DR2Rev, WISE_W1, PS_i, Tycho_V, SDSS_i, WISE_W4, Gaia_G_EDR3, Gaia_RP_EDR3, SkyMapper_r, Gaia_BP_MAWb, SDSS_z, Tycho_B, Bessell_B, DECam_i, Gaia_G_DR2Rev, DECam_Y, 2MASS_J, Kepler_Kp, DECam_u, GALEX_FUV, GALEX_NUV, PS_open, SDSS_u, DECam_r, SkyMapper_g, SkyMapper_z, PS_g, Kepler_D51, DECam_g, Bessell_R, DECam_z, Hipparcos_Hp, WISE_W2, TESS, 2MASS_H, WISE_W3, SkyMapper_i, Gaia_G_MAW, Gaia_BP_MAWf`.

> **Gaia special case**: If *all three* Gaia EDR3 bands (`Gaia_G_EDR3`, `Gaia_BP_EDR3`, `Gaia_RP_EDR3`) are present **and** have uncertainties, a color-dependent extinction law is used for those bands; otherwise fixed extinction coefficients from the table below are used.

### Optional spectroscopy and priors

* **Effective temperature:** `Teff` with `Teff_ERR` (or any allowed suffix). A minimum systematic floor of **100 K** is added internally.
* **Surface gravity:** `logg` with `logg_ERR`. A minimum **0.1 dex** floor is added.
* **Metallicity prior:** Either `FEH_CAL` or `FEH` with matching error column. A minimum **0.1 dex** floor is added. If neither is present, \[Fe/H] is free within global bounds.
* **Parallax prior:** `PARALLAX` with `PARALLAX_ERR`. Optional zero-point correction in `PARALLAX_ZPC` is **added** to `PARALLAX` before use. Only positive effective parallaxes are used; otherwise the sampler falls back to a wide distance prior.
* **Reddening prior:** `EBV` with optional `EBV_ERR`. If `EBV` is present without an error, an uncertainty of **0.3** mag is assumed.

### Photometry error floor

For each band’s magnitude uncertainty, a **0.1 mag** systematic floor is added internally.

---

## How the fit works (parameters & priors)

The sampler explores a 5D parameter vector:

* **Mass** \[M☉], prior \~ Salpeter (α = 2.35) within `[0.1, 5.0]` and model’s maximum mass for the chosen age/metallicity
* **log10(Age/yr)** uniform in `[5.0, 10.1139]`
* **\[Fe/H]**: truncated normal around the provided prior (if any), else uniform in `[-2.9, 0.9]`
* **Distance** \[pc]: from **parallax prior** if provided; otherwise log-uniform in `[1e3, 2e5]`
* **E(B−V)**: uniform in `[EBV_MIN, EBV_MAX]`, default `[0.0, 1.5]`, or tightened to ±5σ around the `EBV` prior (clipped to `[0, 1.5]`).

**Likelihood** combines photometry (with extinction and distance modulus) and, when available, spectroscopic terms for `Teff` and `logg`.

---

## Behavior in common scenarios

* **Photometry-only mode**: Triggered when *both* `Teff` **and** `logg` are **not** simultaneously available with errors. The fit still runs using photometry (and any priors) and reports posteriors. This is the default for most catalog-only use cases.

* **Missing parallax**: No `PARALLAX` prior → distance prior becomes wide log-uniform `[1e3, 2e5]` pc. If `PARALLAX` is present but ≤0 after adding `PARALLAX_ZPC`, it is ignored and the wide prior is used.

* **Missing EBV**: No `EBV` column → E(B−V) prior is uniform over the global `ebv_range` (default `[0.0, 1.5]`). If `EBV` exists **without** an error, we assume `EBV_ERR = 0.3` mag and set the prior to ±5σ around that value (clipped to `[0, 1.5]`).

* **Photometry naming mismatches**: Columns must exactly match supported *minimint* band names. If your catalog uses different names, rename columns before running (see snippet below).

* **Gaia extinction**: If Gaia triplet (`G,BP,RP`) is available with errors, a color-dependent law is iteratively applied; otherwise scalar `A_λ/E(B−V)` coefficients are used.

---

## Outputs

For each star, the following columns are appended:

* `mass_p16`, `mass_p50`, `mass_p84`, `mass_multimodal`
* `logage_p16`, `logage_p50`, `logage_p84`, `logage_multimodal`
* `feh_p16`, `feh_p50`, `feh_p84`, `feh_multimodal`
* `dist_pc_p16`, `dist_pc_p50`, `dist_pc_p84`, `dist_multimodal`
* `ebv_p16`, `ebv_p50`, `ebv_p84`, `ebv_multimodal`
* `lnZ`, `lnZ_err` (Bayesian evidence and its uncertainty from dynesty)

**Multimodality flags** are binary and computed by KDE peak counting on the 1D posterior for each parameter (requires ≥200 samples; peak prominence ≳ 5% of max). `1` indicates multiple significant modes.

If `debug=True`, per-star artifacts are saved under `output_path/<source_id>/`:

* `runplot.png`, `traceplot.png`, `corner.png`
* `isochrone_check.png` (only if Gaia G and Teff are present)
* `summary.json` containing mode, bands used, p50 values, and lnZ

---

## Extinction coefficients (A\_λ / E(B−V))

When not using the Gaia color-dependent law, these coefficients are applied:

* Gaia EDR3: `G=0.83627×3.1`, `BP=1.08337×3.1`, `RP=0.63439×3.1`
* DECam: `g=3.451`, `r=2.646`, `i=2.103`, `z=1.575`, `Y=1.515`
* SDSS: `u=4.871`, `g=3.560`, `r=2.681`, `i=2.400`, `z=1.899`
* 2MASS: `J=0.987`, `H=0.531`, `Ks=0.164`
* SkyMapper: `u=5.017`, `v=4.750`, `g=3.651`, `r=3.032`, `i=2.325`, `z=1.790`
* WISE: `W1=W2=W3=W4=0.0`
* Johnson–Cousins/Bessell: `U=5.47`, `B=4.32`, `V=3.31`, `R=2.68`, `I=1.85`

---

## Tuning & performance

* **`nlive`**: Initial live points for dynesty’s dynamic sampler (`nlive_init`). Larger values better resolve multi-modality and tails, at higher cost. Typical: 2000–8000.
* **`dlogz`**: Evidence tolerance (`dlogz_init`). Smaller is more thorough.
* **`processes`**: Number of processes in an internal `multiprocessing.Pool`. One nested run is executed per star, but dynesty also parallelizes within a run.
* **`random_seed`**: Seeds the RNG for reproducibility across stars.
* **`ebv_range`**: Global E(B−V) bounds when no star-specific EBV prior is available.

Checkpoints (`*_checkpoint.h5`) are written only if `debug=True`.

---

## End-to-end minimal schema examples

### Photometry-only (Gaia + DECam)

Required (any ≥3 bands with errors):

```
source_id
Gaia_G_EDR3, Gaia_G_EDR3_ERR
Gaia_BP_EDR3, Gaia_BP_EDR3_ERR
Gaia_RP_EDR3, Gaia_RP_EDR3_ERR
DECam_g, DECam_g_ERR
```

Optional priors:

```
EBV, EBV_ERR
PARALLAX, PARALLAX_ERR, PARALLAX_ZPC
FEH (or FEH_CAL), FEH_ERR
```

### Photometry + spectroscopy

Add:

```
Teff, Teff_ERR
logg, logg_ERR
```

### Renaming helper (example)

```python
# Suppose your catalog has G_BP, G_RP, G (Gaia EDR3). Rename to minimint names:
for old, new in {
    'phot_g_mean_mag': 'Gaia_G_EDR3',
    'phot_bp_mean_mag': 'Gaia_BP_EDR3',
    'phot_rp_mean_mag': 'Gaia_RP_EDR3',
}.items():
    if old in tab.colnames and new not in tab.colnames:
        tab.rename_column(old, new)
# Likewise ensure error columns exist with one of the accepted suffixes.
```

---

## API reference

```python
fit_stars_with_minimint(
    table: astropy.table.Table,
    output_path: str,
    nlive: int = 3000,
    dlogz: float = 0.01,
    processes: int = 4,
    debug: bool = False,
    random_seed: int = 42,
    ebv_range: tuple[float, float] = (0.0, 1.5),
) -> Table
```

**Returns** the input `Table` with posterior and diagnostic columns appended, and writes an updated table to `output_path`.

---

## Installation & requirements

* Python ≥ 3.9 recommended
* Dependencies: `numpy`, `scipy`, `matplotlib` (headless OK), `astropy`, `dynesty`, `minimint`, `corner` (optional, for debug corner plots)

Install with pip (example):

```bash
pip install numpy scipy matplotlib astropy dynesty corner
# minimint is a separate package; install from pip or source as appropriate
pip install minimint
```

---

## Notes & caveats

* Photometric/spectroscopic errors get minimum floors (0.1 mag, 100 K, 0.1 dex) to account for systematics.
* Distance prior reverts to broad log-uniform if parallax is unusable (non-positive after ZP correction).
* Evidence (`lnZ`) can be compared between phot-only vs spec+phot runs for model comparison, but be mindful of different likelihood terms.
* Multimodality flags are heuristic; inspect posteriors/plots for complex cases.

---

## License

Do whatever you want, fully open source

