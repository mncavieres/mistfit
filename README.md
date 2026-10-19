# mistfit

Bayesian stellar parameters and distances from photometry, by fitting
[MIST](https://waps.cfa.harvard.edu/MIST/) isochrones with the
[`minimint`](https://github.com/segasai/minimint) interpolator and the
[`dynesty`](https://github.com/joshspeagle/dynesty) nested sampler.

Give it an `astropy` table of stars — at least three photometric bands with
uncertainties, plus whatever else you have (Gaia parallax, E(B−V),
spectroscopic Teff / logg / [Fe/H]) — and it returns the same table with
posterior percentiles for **mass, age, [Fe/H], distance and E(B−V)**, the
Bayesian evidence, and multimodality flags for every star.

```
5 free parameters per star:  mass · log10(age) · [Fe/H] · distance · E(B−V)
likelihood:                  photometry (reddened, at distance) + Teff/logg (optional)
                             + parallax (optional, see PARALLAX_MODE)
```

---

## Installation

```bash
pip install mistfit
```

Python ≥ 3.9. This pulls in `numpy`, `scipy`, `astropy`, `matplotlib`,
`tqdm`, `dynesty` and `minimint`. For corner plots in debug mode, also
`pip install "mistfit[plot]"`.

### One-time step: download the MIST grid

`minimint` needs the MIST isochrones and bolometric-correction tables on disk,
and `pip` does not fetch them. Do this once (it downloads a few GB and takes a
while):

```python
import minimint
minimint.download_and_prepare()   # default filter systems: Gaia, 2MASS, Bessell, DECam, SDSS, SkyMapper, WISE, ...
```

Only the filter systems you prepare can be fitted. The VISTA bands are not in
the default list; if you have VHS / VIKING / VVV photometry, add them:

```python
minimint.download_and_prepare(filters=["UBVRIplus", "DECam", "SDSSugriz",
                                       "SkyMapper", "WISE", "VISTA"])
```

Set `MINIMINT_DATA_PATH` to choose where the grid lives (useful on a cluster,
where it can be prepared once and shared).

---

## Quick start

### Minimal working example

Two stars with synthetic photometry generated from MIST itself, so you know
the right answer: a main-sequence dwarf at **1.5 kpc** with Gaia + 2MASS and a
good parallax, and a metal-poor red giant at **20 kpc** whose parallax is
pure noise but which has a spectrum.

```python
# mwe.py
import numpy as np
from astropy.table import Table
from mistfit import fit_stars_with_minimint

if __name__ == "__main__":   # required: the fit runs in a multiprocessing pool
    tab = Table({
        "source_id":    [1, 2],
        "Gaia_G_EDR3":  [15.371, 14.787], "Gaia_G_EDR3_ERR":  [0.003, 0.003],
        "Gaia_BP_EDR3": [15.673, 15.478], "Gaia_BP_EDR3_ERR": [0.005, 0.005],
        "Gaia_RP_EDR3": [14.916, 13.983], "Gaia_RP_EDR3_ERR": [0.005, 0.005],
        "2MASS_J":      [14.402, 12.893], "2MASS_J_ERR":      [0.03, 0.02],
        "2MASS_H":      [14.076, 12.111], "2MASS_H_ERR":      [0.04, 0.02],
        "2MASS_Ks":     [14.032, 11.991], "2MASS_Ks_ERR":     [0.05, 0.02],
        "PARALLAX":     [0.66, 0.03],     "PARALLAX_ERR":     [0.02, 0.05],   # mas
        "EBV":          [0.05, 0.03],     "EBV_ERR":          [0.01, 0.01],   # e.g. from a dust map
        # spectroscopy for star 2 only; NaN means "not measured"
        "Teff": [np.nan, 4190.0], "Teff_ERR": [np.nan, 100.0],
        "logg": [np.nan, 1.1],    "logg_ERR": [np.nan, 0.2],
        "FEH":  [np.nan, -1.0],   "FEH_ERR":  [np.nan, 0.1],
    })

    out = fit_stars_with_minimint(tab, output_path="mistfit_out",
                                  nlive=500, dlogz=0.5, processes=4)

    out["source_id", "fit_mode", "dist_pc_p16", "dist_pc_p50", "dist_pc_p84",
        "feh_p50", "ebv_p50"].pprint()
```

Run it with the giant-friendly sampling coordinate (see
[Configuration](#configuration)):

```bash
SAMPLE_COORD=mass_age_weighted python mwe.py
```

It takes a few minutes on a laptop and prints something like

```
source_id fit_mode dist_pc_p16 dist_pc_p50 dist_pc_p84 feh_p50 ebv_p50
--------- -------- ----------- ----------- ----------- ------- -------
        1     phot      1455.2      1495.8      1544.6  -0.173   0.045
        2     spec     18726.5     21210.6     25462.5  -0.856   0.075
```

Both distances recover the truth (1500 pc and 20 000 pc). The dwarf's age and
[Fe/H] are left essentially unconstrained: broadband photometry alone cannot
pin them on the main sequence. Exact numbers vary slightly between machines
and package versions. `nlive=500, dlogz=0.5` is for a quick demonstration;
use the defaults (`nlive=3000, dlogz=0.01`) or more for science.

> **The `if __name__ == "__main__":` guard is not optional.** The fit uses a
> `multiprocessing.Pool`, and on macOS (and Windows) worker processes re-import
> your script. Without the guard every worker starts its own fit and the run
> hangs with no output and ~0% CPU.

> **Why the E(B−V) and [Fe/H] priors?** Without them, the giant above comes out
> at ~10 kpc with E(B−V) ≈ 0.24: reddening, metallicity and distance are
> partly degenerate in broadband colours. Supply any prior you have.

### A real catalogue

```python
from astropy.table import Table
from mistfit import fit_stars_with_minimint

if __name__ == "__main__":
    tab = Table.read("my_stars.fits")

    # rename your columns to minimint band names (see "Input table format")
    for old, new in {"phot_g_mean_mag": "Gaia_G_EDR3",
                     "phot_bp_mean_mag": "Gaia_BP_EDR3",
                     "phot_rp_mean_mag": "Gaia_RP_EDR3",
                     "parallax": "PARALLAX", "parallax_error": "PARALLAX_ERR"}.items():
        if old in tab.colnames:
            tab.rename_column(old, new)
    # ...and provide an error column for every band, e.g. Gaia_G_EDR3_ERR

    out = fit_stars_with_minimint(tab, output_path="fits/", processes=8)
    # results are also written to fits/fit_results.fits
```

### Command line

The same fit, from a file, with no Python:

```bash
mistfit my_stars.fits --outdir fits/ --nlive 3000 --dlogz 0.01 --procs 8
```

Any table format `astropy` can read works (FITS, ECSV, …). The output is
`fits/fit_results.fits` (or `.ecsv` if the column types cannot be written to
FITS). Add `--debug` for per-star diagnostic plots. Settings from
[Configuration](#configuration) go in the environment:

```bash
SAMPLE_COORD=mass_age_weighted PARALLAX_MODE=likelihood DISTANCE_PRIOR=flat \
  mistfit my_stars.fits --outdir fits/ --procs 8
```

`mistfit --help` (or `python -m mistfit.cli --help`) lists the options.

---

## Input table format

**Minimum requirement**: `≥ 3` photometric bands recognized by *minimint* **and** their per-band uncertainties, using the naming conventions below. Each row corresponds to one star. A `source_id` column is recommended (used to name per-star debug folders), but not required.

### Column names and uncertainty suffixes

* **Photometric band columns** use the *minimint* band names (see list below). Examples: `Gaia_G_EDR3`, `DECam_g`, `2MASS_Ks`.
* The **uncertainty column** for a quantity `<NAME>` is discovered by trying these suffixes (first match wins):

  * `<NAME>_ERR`, `<NAME>_ERRMAG`, `<NAME>_SIG`, `<NAME>_E`, `<NAME>_error`

  Examples: `DECam_g_ERR`, `Gaia_BP_EDR3_ERR`, `Teff_ERR`, `PARALLAX_E`.

### Supported photometric bands

The code auto-detects any of the following present in your table (case-sensitive):
`Bessell_I, PS_z, 2MASS_Ks, SkyMapper_u, Gaia_RP_MAW, Gaia_RP_DR2Rev, Gaia_BP_EDR3, PS_w, PS_r, SDSS_g, SkyMapper_v, Bessell_V, PS_y, SDSS_r, Bessell_U, Gaia_BP_DR2Rev, WISE_W1, PS_i, Tycho_V, SDSS_i, WISE_W4, Gaia_G_EDR3, Gaia_RP_EDR3, SkyMapper_r, Gaia_BP_MAWb, SDSS_z, Tycho_B, Bessell_B, DECam_i, Gaia_G_DR2Rev, DECam_Y, 2MASS_J, Kepler_Kp, DECam_u, GALEX_FUV, GALEX_NUV, PS_open, SDSS_u, DECam_r, SkyMapper_g, SkyMapper_z, PS_g, Kepler_D51, DECam_g, Bessell_R, DECam_z, Hipparcos_Hp, WISE_W2, TESS, 2MASS_H, WISE_W3, SkyMapper_i, Gaia_G_MAW, Gaia_BP_MAWf`.

> **Bands with no extinction coefficient are excluded.** 23 of the names above
> have no entry in `EXT_COEFF` and so cannot be de-reddened: every `PS_*` band,
> `DECam_u`, `TESS`, `Tycho_B`, `Tycho_V`, `Hipparcos_Hp`, `GALEX_FUV`,
> `GALEX_NUV`, `Kepler_Kp`, `Kepler_D51` and the Gaia `DR2Rev`/`MAW` variants.
> Before the 2026-09 release these were accepted, reported as used, counted
> towards the three-band minimum, and then silently ignored by the likelihood.
> They are now excluded with a `MistfitWarning`. Add a coefficient to
> `EXT_COEFF` to bring one back. `VISTA_Y/J/H/Ks` **are** supported, with VHS
> coefficients — do not pass VISTA photometry under 2MASS column names.

> **Gaia special case**: If *all three* Gaia EDR3 bands (`Gaia_G_EDR3`, `Gaia_BP_EDR3`, `Gaia_RP_EDR3`) are present **and** have uncertainties, a color-dependent extinction law is used for those bands; otherwise fixed extinction coefficients from the table below are used.

### Optional spectroscopy and priors

* **Effective temperature:** `Teff` with `Teff_ERR` (or any allowed suffix). **100 K** is *added* to the quoted error (`TEFF_ERR_ADDITIVE`); an optional `TEFF_ERR_FLOOR` then imposes a minimum on the result.
* **Surface gravity:** `logg` with `logg_ERR`. **0.1 dex** added (`LOGG_ERR_ADDITIVE`), optional `LOGG_ERR_FLOOR`.
* **Metallicity prior:** Either `FEH_CAL` or `FEH` with matching error column. **0.1 dex** added (`FEH_ERR_ADDITIVE`), optional `FEH_ERR_FLOOR`. If neither is present, \[Fe/H] is free within global bounds.
* **Parallax prior:** `PARALLAX` with `PARALLAX_ERR`. Optional zero-point correction in `PARALLAX_ZPC` is **added** to `PARALLAX` before use. Since 2026-09 the measurement is used **whatever its sign** (`PARALLAX_MODE`); the old behaviour of discarding non-positive parallaxes is `PARALLAX_MODE=published`.
* **Reddening prior:** `EBV` with optional `EBV_ERR`. If `EBV` is present without an error, an uncertainty of **0.1** mag is assumed.

### Photometry error budget

**0.1 mag** is *added* to each band’s quoted uncertainty (`PHOT_ERR_ADDITIVE`).
An optional `PHOT_ERR_FLOOR` then imposes a minimum on the result. Note these were always
additive despite being described as floors in earlier versions of this file;
both are now available and independently configurable.

---

---

## Outputs

For each star, the following columns are appended:

* `mass_p16`, `mass_p50`, `mass_p84`, `mass_multimodal`
* `logage_p16`, `logage_p50`, `logage_p84`, `logage_multimodal`
* `feh_p16`, `feh_p50`, `feh_p84`, `feh_multimodal`
* `dist_pc_p16`, `dist_pc_p50`, `dist_pc_p84`, `dist_multimodal`
* `ebv_p16`, `ebv_p50`, `ebv_p84`, `ebv_multimodal`
* `lnZ`, `lnZ_err` (Bayesian evidence and its uncertainty from dynesty)
* `fit_mode`: `spec` if both `Teff` and `logg` entered the likelihood, else `phot`

**Multimodality flags** are binary and computed by KDE peak counting on the 1D posterior for each parameter (requires ≥200 samples; peak prominence ≳ 5% of max). `1` indicates multiple significant modes.

If `debug=True`, per-star artifacts are saved under `output_path/<source_id>/`:

* `runplot.png`, `traceplot.png`, `corner.png`
* `isochrone_check.png` (only if Gaia G and Teff are present)
* `summary.json` containing mode, bands used, p50 values, and lnZ

---

## Configuration

Modelling choices are set through **environment variables**, read when
`mistfit` is imported. That makes it easy to run the same script under a job
scheduler with one thing changed. From Python, set them before the import (the
worker processes inherit the environment):

```python
import os
os.environ["SAMPLE_COORD"] = "mass_age_weighted"
os.environ["LOGAGE_MIN"], os.environ["LOGAGE_MAX"] = "9.6", "10.13"   # 4-13.5 Gyr
from mistfit import fit_stars_with_minimint
```

| variable | default | what it does |
|---|---|---|
| `SAMPLE_COORD` | `mass_age` | sampling coordinates: `mass_age`, `mass_eep` or `mass_age_weighted`. **Use `mass_age_weighted` for giants.** |
| `PARALLAX_MODE` | `prior` | how the parallax enters: `prior` (sets the distance prior), `likelihood` (Gaussian term on 1/d), or `published` (old behaviour; discards ϖ ≤ 0) |
| `DISTANCE_PRIOR` | `loguniform` | prior on distance when the parallax is not acting as one: `loguniform`, `flat` or `volume` (∝ d²) |
| `LOGAGE_MIN`, `LOGAGE_MAX` | `5.0`, `10.114` | age window, log10(yr). The default lower limit admits pre-main-sequence solutions; for old populations narrow it. |
| `EEP_MIN`, `EEP_MAX` | `200`, `808` | EEP window under `SAMPLE_COORD=mass_eep` |
| `PHOT_ERR_ADDITIVE` | `0.1` mag | added to every photometric error |
| `TEFF_ERR_ADDITIVE` | `100` K | added to the Teff error |
| `LOGG_ERR_ADDITIVE` | `0.1` dex | added to the logg error |
| `FEH_ERR_ADDITIVE` | `0.1` dex | added to the [Fe/H] error |
| `*_ERR_FLOOR` | none | minimum error after the addition (`PHOT_`, `TEFF_`, `LOGG_`, `FEH_`) |
| `RNG_PER_STAR` | `1` | seed each star from (`random_seed`, `source_id`), so a star's result does not depend on the rest of the table |

Function arguments (`nlive`, `dlogz`, `processes`, `random_seed`, `ebv_range`,
`debug`) are described under [API reference](#api-reference).

### Choosing settings for your science case

There is no single right answer here, which is why these are settings rather
than a decision baked into the code. What follows is guidance, not a
recommendation — the right choice depends on how informative your parallaxes
are and where your targets actually sit.

#### The one question that decides most of it

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
| reproducing a pre-2026-09 run | `PARALLAX_MODE=published RNG_PER_STAR=0 SAMPLE_COORD=mass_age` | see CHANGELOG; it will not be bit-identical |

#### The other question: are your stars evolved?

**`SAMPLE_COORD` decides how well the sampler can resolve the giant branch.**
Under the default `mass_age`, at a fixed old age the whole RGB lies in a narrow
sliver of mass just below `getMaxMass(logAge, feh)`, and the likelihood cuts it
off there. Two ways round that:

* `mass_eep` samples EEP instead of age, so the boundary disappears: the
  constraint is simply whether (mass, EEP, [Fe/H]) exists on the MIST grid. The
  price is the prior — uniform in evolutionary stage, not in age.
* `mass_age_weighted` keeps (mass, logAge) and draws the mass from a proposal
  crowded against `getMaxMass`, dividing it back out of the likelihood. Same
  prior, posterior and `lnZ` as `mass_age`; only the sampling changes.

| your sample | suggested | why |
|---|---|---|
| giants (RGB, AGB), or anything where `logg ≲ 2` | `SAMPLE_COORD=mass_age_weighted` | reaches the giant branch and the isochrone's end with `mass_age`'s prior, so you know what the prior is; on 51 K giants it found 10 and 37 nats more summed evidence than `mass_age` (same model) on the two samples, had the fewest stars in >3σ tension with their spectroscopic `logg` (3, against 7 and 5), and used the fewest likelihood calls of the three |
| giants, where you want MINESweeper's convention | `SAMPLE_COORD=mass_eep` | uniform in EEP; read the caveats below before comparing its distances or `lnZ` with the other two |
| dwarfs and subgiants | any | the main sequence is well resolved in all three coordinates; on DESI312 they agree within 0.16σ |
| reproducing an older run | `SAMPLE_COORD=mass_age` | the default, unchanged |

**What uniform-in-EEP does to the prior.** It is not a reparameterisation of
`mass_age`; it is a different prior, in two ways:

* inside its support it weights each phase by how fast the star is changing
  (dEEP/dlogAge), not by how long it stays there — the RGB at `logg = 1` gets
  ~550× the main sequence's weight relative to a uniform age prior, and the
  weight changes by a factor 2.6 across `logg = 2.2 ± 0.3`;
* its support ends at `EEP_MAX = 808`, the start of the TP-AGB, so TP-AGB and
  post-AGB solutions are excluded outright. On 9 of 51 K giants most of the
  uniform-age posterior sat there, and `mass_eep` put those stars 11–36%
  nearer. Raising `EEP_MAX` swaps one problem for another: the TP-AGB is 601
  EEPs long.

The ~85 nats by which `mass_eep` beat `mass_age` on 18 MagE giants are mostly a
Bayes factor between those two priors (74.6 nats), not better sampling (10.3).
Details in CHANGELOG 2026-09c.

The output columns are identical in all three — the EEP is converted back to
logAge before anything is summarised. Only a raw dynesty **checkpoint** under
`mass_eep` holds the EEP in column 1.

#### The trade-off between `prior` and `likelihood`

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

#### Error floors

Use `*_ERR_FLOOR` when your pipeline reports uncertainties that are implausibly
small rather than merely optimistic — spectroscopic `logg` quoted to 0.0004 dex
will otherwise pin the fit wherever it landed. Use `*_ERR_ADDITIVE` (the
default) when the quoted errors are believable and you want to fold in a
systematic. They compose: `max(err + additive, floor)`.

Note the interaction: a loose gravity constraint and a discarded parallax
together are what let fits collapse onto pre-main-sequence solutions. Either
one alone was enough to prevent it in testing, so if you must use a wide `logg`
floor, keep the parallax.

#### Diagnostics worth running whatever you choose

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

---

## How the fit works (parameters & priors)

The sampler explores a 5D parameter vector:

* **Mass** \[M☉], prior \~ Salpeter (α = 2.35) within `[0.1, 100]` (`M_MIN`, `M_MAX`) and, under `SAMPLE_COORD=mass_age`, below the model’s maximum mass for the chosen age/metallicity
* **log10(Age/yr)** uniform in `[LOGAGE_MIN, LOGAGE_MAX]`, default `[5.0, 10.1139]` (13 Gyr), both settable from the environment. Under `SAMPLE_COORD=mass_eep` logAge is derived from the EEP instead, and points outside the window are rejected.
* **\[Fe/H]**: truncated normal around the provided prior (if any), else uniform in `[-2.9, 0.9]`
* **Distance** \[pc], within `[10, 2e5]` (`DIST_MIN`, `DIST_MAX`): from the **parallax** if one is given — `1000/ϖ`, with ϖ drawn from the measurement's Gaussian truncated at zero — otherwise from `DISTANCE_PRIOR`, log-uniform by default. Under `PARALLAX_MODE=likelihood` the distance always comes from `DISTANCE_PRIOR` and the parallax enters the likelihood instead.
* **E(B−V)**: uniform over `ebv_range`, default `[EBV_MIN_DEFAULT, EBV_MAX_DEFAULT]` = `[0.0, 1.5]`; for a star with an `EBV` value, over ±5σ around it, clipped to `[EBV_MIN_DEFAULT, EBV_MAX_DEFAULT]`.

**Likelihood** combines photometry (with extinction and distance modulus), spectroscopic terms for `Teff` and `logg` when both are available, and the parallax under `PARALLAX_MODE=likelihood`.

---

## Behavior in common scenarios

* **Photometry-only mode**: Triggered when *both* `Teff` **and** `logg` are **not** simultaneously available with errors. The fit still runs using photometry (and any priors) and reports posteriors. This is the default for most catalog-only use cases.

* **Missing parallax**: No `PARALLAX` prior → the distance prior is whatever `DISTANCE_PRIOR` says, over `[DIST_MIN, DIST_MAX]`. A present but non-positive `PARALLAX` is **used** since 2026-09 (`PARALLAX_MODE`), not ignored — for a distant sample a negative measured parallax is noise on a small positive number, and discarding it is a selective loss.

* **Missing EBV**: No `EBV` column → E(B−V) prior is uniform over the global `ebv_range` (default `[0.0, 1.5]`). If `EBV` exists **without** an error, we assume `EBV_ERR = 0.1` mag and set the prior to ±5σ around that value (clipped to `[EBV_MIN_DEFAULT, EBV_MAX_DEFAULT]`, default `[0, 1.5]`).

* **Photometry naming mismatches**: Columns must exactly match supported *minimint* band names. If your catalog uses different names, rename columns before running (see [A real catalogue](#a-real-catalogue)).

* **Gaia extinction**: If Gaia triplet (`G,BP,RP`) is available with errors, a color-dependent law is iteratively applied; otherwise scalar `A_λ/E(B−V)` coefficients are used.

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
* VISTA (VHS): `Y=1.310`, `J=0.928`, `H=0.592`, `Ks=0.388`

---

## Performance tuning

* **`nlive`**: Initial live points for dynesty’s dynamic sampler (`nlive_init`). Larger values better resolve multi-modality and tails, at higher cost. Typical: 2000–8000.
* **`dlogz`**: Evidence tolerance (`dlogz_init`). Smaller is more thorough.
* **`processes`**: Number of processes in an internal `multiprocessing.Pool`. One nested run is executed per star, but dynesty also parallelizes within a run.
* **`random_seed`**: Seeds the RNG for reproducibility across stars.
* **`ebv_range`**: Global E(B−V) bounds when no star-specific EBV prior is available.

Checkpoints (`*_checkpoint.h5`) are written only if `debug=True`.

---

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
    ebv_range: tuple[float, float] | None = None,  # None: (EBV_MIN_DEFAULT, EBV_MAX_DEFAULT) = (0.0, 1.5)
) -> Table
```

**Returns** the input `Table` with posterior and diagnostic columns appended, and writes an updated table to `output_path`.

---

---

## Notes & caveats

* Photometric/spectroscopic errors are inflated additively by default (0.1 mag, 100 K, 0.1 dex); optional floors are available. See `CHANGELOG.md`.
* The Gaia parallax is used regardless of sign by default, and the distance prior used when it is not acting as one is selectable (`DISTANCE_PRIOR`). The historical `loguniform` default puts **70% of its mass inside 10 kpc** over the default bounds, which matters if your targets are distant.
* Evidence (`lnZ`) is comparable only between fits with the same likelihood — `PARALLAX_MODE=published` against `prior`, say, or one `DISTANCE_PRIOR` against another. A spec+phot and a phot-only fit of the same star are not comparable, and neither is `PARALLAX_MODE=likelihood` against the other modes: each adds likelihood terms, so its `lnZ` is the evidence for different data.
* Multimodality flags are heuristic; inspect posteriors/plots for complex cases.

---

## Upgrading from a pre-2026-09 version

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
| **`SAMPLE_COORD`** | `mass_age` \| `mass_eep` \| `mass_age_weighted`. Sampling in EEP rather than age reaches the giant branch, which spans 47.7% of a U(200, 808) EEP prior against 1.8% of a 4–13 Gyr age prior. `mass_age_weighted` reaches it too while keeping `mass_age`'s prior and evidence (2026-09c). |

Configuration is read from the environment, so a scheduler run can change one
setting and nothing else:

```bash
PARALLAX_MODE=prior DISTANCE_PRIOR=volume LOGG_ERR_FLOOR=0.3 python run_nest.py
```

---

## Citing mistfit

If you use `mistfit` in published work, please cite:

> Cavieres, M., Koposov, S. E., Rossi, E. M., Penoyre, Z., & Verberne, S. 2026,
> *Discovery of Galactic-centre-ejected star in DESI DR1*, A&A, 709, A117.
> [ADS](https://ui.adsabs.harvard.edu/abs/2026A%26A...709A.117C/abstract) ·
> [doi:10.1051/0004-6361/202659104](https://doi.org/10.1051/0004-6361/202659104) ·
> [arXiv:2601.19866](https://arxiv.org/abs/2601.19866)

```bibtex
@ARTICLE{Cavieres2026,
       author = {{Cavieres}, Manuel and {Koposov}, Sergey E. and {Rossi}, Elena Maria and
                 {Penoyre}, Zephyr and {Verberne}, Sill},
        title = "{Discovery of Galactic-centre-ejected star in DESI DR1}",
      journal = {\aap},
         year = 2026,
       volume = {709},
          eid = {A117},
          doi = {10.1051/0004-6361/202659104},
       eprint = {2601.19866},
archivePrefix = {arXiv},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2026A%26A...709A.117C}
}
```

`mistfit` is a thin layer over other people's work; please also cite:

* **MIST** — [Dotter 2016, ApJS, 222, 8](https://doi.org/10.3847/0067-0049/222/1/8);
  [Choi et al. 2016, ApJ, 823, 102](https://doi.org/10.3847/0004-637X/823/2/102)
* **dynesty** — [Speagle 2020, MNRAS, 493, 3132](https://doi.org/10.1093/mnras/staa278)
* **minimint** — [Koposov, github.com/segasai/minimint](https://github.com/segasai/minimint)

---

## License

MIT — see [LICENSE](LICENSE).
