# Changelog

## 2026-09 — four bugs fixed, two modelling choices exposed

Three of these change results. Read "What this changes for you" at the bottom
before rerunning anything you have already published.

---

### FIXED — the band list was iterated in hash order

`_bands_for_row` built its list by iterating `MINIMINT_BANDS`, which is a
`set`. Python randomises string hashing per process, so the order of the
returned list varied between runs. That list is handed to
`loglike_phot_theta`, which accumulates

```python
ll += -0.5 * ((o - (mod + dm + A)) / e)**2
```

over it. Floating-point addition is not associative, so two runs of the same
star on the same data got log-likelihoods differing in the last bits — and on
posteriors as multimodal as these, that is enough to flip which mode the
sampler settles in.

Measured on a 33-star sample fitted twice with identical code, identical input
and the same `random_seed`: **every one of the 33 stars differed by more than
1%, 21 by more than 10%, with distance ratios from 0.37 to 4.63** and one
star's `lnZ` moving by 72. The two runs' logged `obs` dictionaries carry the
same numbers in a different order, which is the fingerprint.

Fixed by `for b in sorted(MINIMINT_BANDS)`. Setting `PYTHONHASHSEED` is no
longer necessary, though it remains harmless.

### FIXED — bands with no extinction coefficient were dropped in silence

`loglike_phot_theta` does `EXT_COEFF.get(band)` and `continue`s on a miss. **23
of the 58 names in `MINIMINT_BANDS` have no `EXT_COEFF` entry**, among them
every Pan-STARRS band, `DECam_u`, `TESS`, `Tycho_B`, `Tycho_V`,
`Hipparcos_Hp`, both GALEX bands and the Gaia DR2Rev/MAW variants. Photometry
in those bands was accepted, reported in `summary.json` as used, counted
towards the `len(bands) >= 3` requirement — and then contributed nothing to the
likelihood.

Anyone fitting Pan-STARRS photometry was getting a fit that ignored it.

Now excluded in `_bands_for_row` with a `RuntimeWarning` naming the band, once
per process. The reported band list is therefore the list that was used, and
the three-band minimum counts only bands that do something. Add an `EXT_COEFF`
entry to bring a band back.

### FIXED — one random stream shared across the whole table

`fit_stars_with_minimint` created a single `np.random.default_rng(random_seed)`
and passed it to every star in turn. The stream a star received therefore
depended on how many draws the stars before it had consumed. Results moved when
one star's sampling diverged, when the table was reordered, and when a run was
split across jobs — so two fits of the same star could not be compared.

Now one `Generator` per star, seeded from `SeedSequence([random_seed,
source_id])`, so a star's stream depends only on which star it is. A run split
into one-star jobs gives the same answer as the whole table in one job.

`RNG_PER_STAR=0` restores the old behaviour for reproducing earlier results.

### FIXED — `__version__` was always `0+unknown`

`__init__.py` asked `importlib.metadata` for the version of
`"mist-stellar-fitter"`, a name the package has not had since it was renamed to
`mistfit`, so the lookup always raised and fell through to the fallback. It also
exported `__all__ = ["mistfit"]`, naming the module rather than the function it
provides. Both fixed.

---

### NEW — `PARALLAX_MODE`: the parallax is no longer discarded when negative

The old code did:

```python
if plx > 0:
    obs['parallax'] = (plx, float(perr))
```

A true parallax is always positive, but the *measurement* scatters either side
of zero, and for a distant star the true value is far smaller than the
uncertainty. At 60 kpc the true parallax is 0.017 mas against a Gaia
uncertainty of 0.04–0.15 mas, so **roughly a third of genuinely distant stars
measure negative on noise alone** — 41% at 100 kpc. In one 18-star sample, 13
measured negative and not one was even 2σ from zero.

Discarding those is a discard *conditional on the noise realisation*, and it is
selective: only the stars that scattered towards "far" lose their distance
prior and fall back on the wide default. It is a defect, not a conservative
choice, which is why the corrected behaviour is now the default.

| value | behaviour |
|---|---|
| `"published"` | the old gate. For reproducing pre-2026-09 results. |
| `"prior"` | **(default)** parallax kept whatever its sign, as the distance prior |
| `"likelihood"` | parallax kept, entering the likelihood as a Gaussian observation of 1/d |

Measured on an 18-star sample where the constraint bites (spectroscopic
gravities given a 0.3 dex floor, so the isochrone fit has room to wander):

| | `published` | `prior` | `likelihood` |
|---|---|---|---|
| collapsed onto pre-main-sequence solutions | 11/18 | **2/18** | 5/18 |
| disagree with their own Gaia parallax by >3σ | 8/18 | **0/18** | 4/18 |
| median \|Δϖ\| | 2.09σ | **0.90σ** | 1.16σ |
| Σ lnZ | −748 | **−463** | −728 |

`prior` against `published` is a legitimate Bayes factor — same likelihood,
different prior — and it is e^285 in favour of keeping the measurement.

**On a sample with a tighter gravity constraint the effect is much smaller.** A
33-star sample using the default additive errors showed 1/33 collapses and 0/33
parallax disagreements under *every* mode; the parallax changed individual
distances but fixed nothing, because the gravity prior was already preventing
the failure. Both conditions have to be present.

**Why `likelihood` is not the default, despite being formally cleaner.** It has
no truncation, no clip and no Jacobian, and it is what MINESweeper and SpecDis
do. But it turns a hard geometric constraint into a soft penalty that a
competing likelihood mode can outbid, and it samples less reliably: on the
18-star sample four stars fell into a pre-main-sequence basin whose `lnZ` was
31–87 nats *worse* than the solution `prior` found for the same star. The two
posteriors differ by exactly one factor of 1/d in the effective prior, worth at
most ln(200) ≈ 5 nats across the whole range, so a gap of 87 is the sampler
failing to find the mode rather than the model preferring another one. If you
use `likelihood`, check `DISTANCE_PRIOR` too — see below.

### NEW — `DISTANCE_PRIOR`

The prior on distance used whenever the parallax is not itself acting as one:
always under `PARALLAX_MODE="likelihood"`, for discarded stars under
`"published"`, and for any star with no usable parallax.

Over the default bounds (`DIST_MIN = 10 pc`, `DIST_MAX = 200 kpc`):

| value | p(d) | median | mass below 10 kpc | weight given to 2 kpc vs 50 kpc |
|---|---|---|---|---|
| `"loguniform"` | ∝ 1/d | 1.4 kpc | 69.8% | 25 |
| `"flat"` | ∝ const | 100 kpc | 5.0% | 1 |
| `"volume"` | ∝ d² | 159 kpc | 0.0% | 0.0016 |

Default stays `"loguniform"`, which is what the code has always done — but note
what that means: **70% of the live points start inside 10 kpc, and the prior
median is 1.4 kpc.** For a sample of distant stars that is a region containing
nothing, and a 2 kpc solution is handed 25 times the prior weight of a 50 kpc
one before any data are seen. If your targets are distant, `"volume"` is worth
trying; it is 9.7 nats less generous to the near solution.

These are prior-shape statements, not recommendations: which one is right
depends on your targets, and `"loguniform"` remains the default precisely so
that nothing changes for you unless you choose it.

One trap worth naming: the obvious "physical" prior for a tracer population,
d²ρ(r) with ρ ∝ r^−3.5, goes as d^−1.5 and favours near stars *more* strongly
than log-uniform. That is correct for a volume-complete sample and wrong for a
magnitude-limited sample of distant giants, where the selection function — not
the density profile — is what makes the sample distant. None of the three
encodes a selection function. `"volume"` is the standard uninformative choice
in three dimensions, not a claim about the Galaxy.

### NEW — configurable error budget, and VISTA bands

`*_ERR_ADDITIVE` and `*_ERR_FLOOR` for Teff, logg, [Fe/H] and photometry. The
additive defaults reproduce the historical +100 K / +0.1 dex / +0.1 dex /
+0.1 mag exactly; floors default to `None`. A floor is the right tool when
quoted errors are implausibly small rather than merely optimistic — a pipeline
reporting logg to 0.0004 dex will otherwise pin the fit wherever it landed.

`VISTA_Y/J/H/Ks` added to `MINIMINT_BANDS` with VHS extinction coefficients.
Near-IR photometry from VHS/VIKING/VVV is routinely passed to fitters under
2MASS column names, which silently matches it against 2MASS bolometric
corrections and de-reddens it with 2MASS coefficients: A_Ks/E(B−V) is 0.388 for
VISTA against 0.164 for 2MASS, a factor 2.4.

---

## What this changes for you

Even with every setting left at its default, **results will move**, because two
of the three bug fixes change numbers:

* the band-order fix makes runs reproducible, but the number you now get
  reproducibly is not necessarily the number you got before;
* excluding bands with no extinction coefficient removes photometry that was
  being accepted and ignored — if you fit Pan-STARRS, your band list changes;
* per-star seeding changes every star's random stream.

To reproduce a pre-2026-09 run as closely as the code now allows:

```bash
PARALLAX_MODE=published RNG_PER_STAR=0 python run_nest.py
```

That still will not be bit-identical, because the band ordering it used was
never reproducible in the first place — which is the point of the fix.
