# Nested Sampling MIST fitter

Nested-sampling stellar-parameter fits using **dynesty** + **minimint** (MIST isochrones with extinction).

This package provides a single public function, `fit_stars_with_minimint`, which takes an `astropy.table.Table` containing photometry (and optionally spectroscopy, parallax, and reddening priors), and returns the same table with posterior summaries appended. A simple CLI is also included.

---
## What changed in 2026-09

Four bugs fixed and two modelling choices exposed as settings. Two of the fixes
change results even with every setting left at its default, so read
[CHANGELOG.md](CHANGELOG.md) before rerunning anything you have published.

| | |
|---|---|
| **band ordering** | the band list was built by iterating a `set`, so the order of the terms summed into the log-likelihood followed `PYTHONHASHSEED`. Two runs of the same stars on the same data gave distances differing by factors of up to 4.6. |
| **silent band drops** | 23 of 58 "supported" bands have no extinction coefficient and were accepted, reported as used, then ignored by the likelihood. Every Pan-STARRS band was among them. |
| **shared random stream** | one `Generator` per table, handed to each star in turn, so results depended on row order and on how the run was split across jobs. |
| **`__version__`** | looked up a package name that no longer existed, so it was always `0+unknown`. |
| **`PARALLAX_MODE`** | the parallax was discarded whenever measured ≤ 0 — which, for distant stars, is noise on a small positive number and a *selective* discard. Now used regardless of sign by default. |
| **`DISTANCE_PRIOR`** | `loguniform` \| `flat` \| `volume`. The historical `loguniform` puts 70% of its prior mass inside 10 kpc. |

Configuration is read from the environment, so a scheduler run can change one
setting and nothing else:

```bash
PARALLAX_MODE=prior DISTANCE_PRIOR=volume LOGG_ERR_FLOOR=0.3 python run_nest.py
```

## Choosing settings for your science case

There is no single right answer here, which is why these are settings rather
than a decision baked into the code. What follows is guidance, not a
recommendation — the right choice depends on how informative your parallaxes
are and where your targets actually sit.

### The one question that decides most of it

**How informative is your parallax?** Not whether it is positive — whether
`|ϖ| / σ_ϖ` is large. For a star at 60 kpc the true parallax is 0.017 mas
against a typical Gaia uncertainty of 0.04–0.15 mas, so the measurement says
little beyond "not nearby". For a star at 2 kpc it is 0.5 mas and pins the
distance on its own.

| your sample | suggested | why |
|---|---|---|
| nearby, ϖ S/N ≳ 5 | `PARALLAX_MODE=likelihood`, `DISTANCE_PRIOR` barely matters | the parallax dominates; the likelihood form is the correct one and the prior is not doing the work |
| distant, ϖ S/N ≲ 2, posteriors multimodal | `PARALLAX_MODE=prior` | the parallax shapes the *proposals*, so the sampler cannot spend live points where the parallax forbids. This is the most robust option when the isochrone posterior has competing modes |
| distant, but you want the cleanest formulation | `PARALLAX_MODE=likelihood`, `DISTANCE_PRIOR=volume` | no truncation, no clip, no Jacobian — but check the diagnostics below, because the likelihood form samples less reliably |
| no parallax at all | `DISTANCE_PRIOR` is doing everything | pick it from where your targets are; `loguniform` puts 70% of its mass inside 10 kpc |
| reproducing a pre-2026-09 run | `PARALLAX_MODE=published RNG_PER_STAR=0` | see CHANGELOG; it will not be bit-identical |

### The trade-off between `prior` and `likelihood`

Both use the parallax regardless of sign. They differ in *where* it enters, and
the difference is not only formal:

* `prior` makes it a hard geometric constraint. The sampler cannot propose a
  distance the parallax excludes. Robust, but it carries an implicit 1/d²
  Jacobian, so for a star whose parallax is genuinely uninformative it can still
  pull the distance inward against the photometry's preference.
* `likelihood` makes it a soft penalty. Formally cleaner — it is what
  MINESweeper and SpecDis do — but a competing likelihood mode can outbid the
  penalty. On one 18-star test set, four stars fell into a pre-main-sequence
  solution whose `lnZ` was 31–87 nats *worse* than what `prior` found for the
  same star.

If you use `likelihood`, pair it with a distance prior that does not seed the
sampler in the wrong place, and check the diagnostics below.

### Error floors

Use `*_ERR_FLOOR` when your pipeline reports uncertainties that are implausibly
small rather than merely optimistic — spectroscopic `logg` quoted to 0.0004 dex
will otherwise pin the fit wherever it landed. Use `*_ERR_ADDITIVE` (the
default) when the quoted errors are believable and you want to fold in a
systematic. They compose: `max(err + additive, floor)`.

Note the interaction: a loose gravity constraint and a discarded parallax
together are what let fits collapse onto pre-main-sequence solutions. Either
one alone was enough to prevent it in testing, so if you must use a wide `logg`
floor, keep the parallax.

### Diagnostics worth running whatever you choose

None of these settings can be validated from the fitted distance alone. Check,
per star:

* how much posterior mass sits below `logage` 6 — a halo giant fitted as a
  pre-main-sequence star is the classic failure;
* how much sits on `DIST_MAX`, where `PARALLAX_MODE=prior` clips;
* the winning model's `logg` against the spectroscopic value you supplied, in σ;
* the parallax implied by the fitted distance against the measured one, in σ,
  *whether or not the run used it*.

A fit that disagrees with its own inputs by many σ is an artifact regardless of
how tight its error bars look.

## Citation

If this tool is used in a publication, please cite **dynesty**, **minimint**, and Cavieres (in prep).


## Quick start

---

## Installation & requirements

* Python ≥ 3.9 recommended
* Dependencies: `numpy`, `scipy`, `matplotlib` (headless OK), `astropy`, `dynesty`, `minimint`, `corner` (optional, for debug corner plots)

Install with pip (example):

```bash
pip install mistfit
```


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
# run a nested fit on an input catalog (FITS/ECSV supported by Astropy)
mistfit /path/to/input_table.fits \
  --outdir /path/to/out \
  --nlive 5000 \
  --dlogz 0.01 \
  --procs 8 \
  --debug

```

The outputs are written to `--outdir`:
--outdir will be created if it doesn’t exist; results are saved as fit_results.fits (or .ecsv fallback).

* Updated table with appended posterior columns: `fit_results.fits` (or `.ecsv` if FITS types conflict)
* If `--debug`, per-star folders with `runplot.png`, `traceplot.png`, `corner.png`, optional `isochrone_check.png`, and `summary.json`.

Helpful commands:

```bash
mistfit --help
```

module form (works even without the console script on PATH)

```bash
python -m mistfit.cli --help
```

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

> **Bands with no extinction coefficient are excluded.** 23 of the names above
> have no entry in `EXT_COEFF` and so cannot be de-reddened: every `PS_*` band,
> `DECam_u`, `TESS`, `Tycho_B`, `Tycho_V`, `Hipparcos_Hp`, `GALEX_FUV`,
> `GALEX_NUV`, `Kepler_Kp`, `Kepler_D51` and the Gaia `DR2Rev`/`MAW` variants.
> Before the 2026-09 release these were accepted, reported as used, counted
> towards the three-band minimum, and then silently ignored by the likelihood.
> They are now excluded with a `RuntimeWarning`. Add a coefficient to
> `EXT_COEFF` to bring one back. `VISTA_Y/J/H/Ks` **are** supported, with VHS
> coefficients — do not pass VISTA photometry under 2MASS column names.

> **Gaia special case**: If *all three* Gaia EDR3 bands (`Gaia_G_EDR3`, `Gaia_BP_EDR3`, `Gaia_RP_EDR3`) are present **and** have uncertainties, a color-dependent extinction law is used for those bands; otherwise fixed extinction coefficients from the table below are used.

### Optional spectroscopy and priors

* **Effective temperature:** `Teff` with `Teff_ERR` (or any allowed suffix). **100 K** is *added* to the quoted error (`TEFF_ERR_ADDITIVE`); an optional `TEFF_ERR_FLOOR` imposes a minimum instead.
* **Surface gravity:** `logg` with `logg_ERR`. **0.1 dex** added (`LOGG_ERR_ADDITIVE`), optional `LOGG_ERR_FLOOR`.
* **Metallicity prior:** Either `FEH_CAL` or `FEH` with matching error column. **0.1 dex** added (`FEH_ERR_ADDITIVE`), optional `FEH_ERR_FLOOR`. If neither is present, \[Fe/H] is free within global bounds.
* **Parallax prior:** `PARALLAX` with `PARALLAX_ERR`. Optional zero-point correction in `PARALLAX_ZPC` is **added** to `PARALLAX` before use. Since 2026-09 the measurement is used **whatever its sign** (`PARALLAX_MODE`); the old behaviour of discarding non-positive parallaxes is `PARALLAX_MODE=published`.
* **Reddening prior:** `EBV` with optional `EBV_ERR`. If `EBV` is present without an error, an uncertainty of **0.3** mag is assumed.

### Photometry error budget

**0.1 mag** is *added* to each band’s quoted uncertainty (`PHOT_ERR_ADDITIVE`).
An optional `PHOT_ERR_FLOOR` imposes a minimum instead. Note these were always
additive despite being described as floors in earlier versions of this file;
both are now available and independently configurable.

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

* **Missing parallax**: No `PARALLAX` prior → the distance prior is whatever `DISTANCE_PRIOR` says, over `[DIST_MIN, DIST_MAX]`. A present but non-positive `PARALLAX` is **used** since 2026-09 (`PARALLAX_MODE`), not ignored — for a distant sample a negative measured parallax is noise on a small positive number, and discarding it is a selective loss.

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

## Notes & caveats

* Photometric/spectroscopic errors are inflated additively by default (0.1 mag, 100 K, 0.1 dex); optional floors are available. See `CHANGELOG.md`.
* The Gaia parallax is used regardless of sign by default, and the distance prior used when it is not acting as one is selectable (`DISTANCE_PRIOR`). The historical `loguniform` default puts **70% of its mass inside 10 kpc** over the default bounds, which matters if your targets are distant.
* Evidence (`lnZ`) can be compared between phot-only vs spec+phot runs for model comparison, but be mindful of different likelihood terms.
* Multimodality flags are heuristic; inspect posteriors/plots for complex cases.

---

## License

Do whatever you want, fully open source

