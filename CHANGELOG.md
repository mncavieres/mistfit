# Changelog

## 2026-09c — a weighted mass proposal for the giant branch

### NEW — `SAMPLE_COORD=mass_age_weighted`

A third sampling coordinate, adapted from `MassMapper` in the rvspecfit-based
distance code. It samples (mass, logAge) exactly as `mass_age` does, but draws
the mass from a proposal that crowds points against the isochrone's own maximum
mass, then divides that proposal back out of the likelihood:

    logL  +=  ln S(m) − ln q(m | logAge, [Fe/H])

where S is the Salpeter density `mass_age` draws from and q is a broken power
law with its edge at `getMaxMass(logAge, [Fe/H])`: flat, then rising as
(1 − m/m_max)^−1.545 over the last 1.4% of the mass range, then flat again
over the last 0.03%. The break points and slope are `MassMapper`'s,
unchanged.

**The model is identical to `mass_age`.** Same prior, same likelihood, same
posterior, and — unlike `mass_eep` — the same evidence, so `lnZ` is directly
comparable between the two. Only where the sampler puts its points changes.
That makes this the way to resolve the giant branch *without* swapping the
uniform age prior for a uniform evolutionary-stage prior, which is the trade
`mass_eep` makes.

What the proposal buys, in the setup behind the 12.8% figure in 2026-09b
([Fe/H] ~ N(−1.2, 0.3), log age over 4–13 Gyr, 4000 prior draws each):

| | on the grid | of those, `logg < 2` | draws at `logg < 2` | lowest `logg` |
|---|---|---|---|---|
| `mass_age` | 3821/4000 | 0.0% | 0 | 2.51 |
| `mass_eep` | 223/4000 | 13.5% | 30 | 0.58 |
| `mass_age_weighted` | **4000/4000** | 12.2% | **490** | **−0.24** |

`mass_eep` reaches the upper giant branch as often per on-grid draw, but most of
its draws fall off the grid; the weighted proposal puts 16× more draws there in
absolute terms. Per isochrone, 25% of proposal draws land in the last 1.4% of
the mass range, where the post-main-sequence sits (0.1% of Salpeter draws),
and 9.6% in the last 0.03%, where the TP-AGB and post-AGB sit (0.002% of
Salpeter draws, a factor of ~4000).

Checks of the weighting itself: `q` is an exact inverse CDF (round trip to
1e−14), its density integrates to one, and `q · exp(weight)` reproduces the
Salpeter density to 1e−10 pointwise and its on-isochrone prior volume to 1e−8.
A 400 000-draw Monte Carlo recovers the giant-branch tail's prior mass to
within its 0.9% sampling error. Tests in `test/test_config.py` pin all of this.

### Measured: a dwarf, and 51 K giants

**DESI312** (a 1 M☉, 10 Gyr main-sequence/turn-off star; parallax in the
likelihood, flat distance prior, nlive 2000). All three coordinates agree
within 0.16σ. `mass_age` and `mass_age_weighted` give lnZ −17.945 and −17.940
— the same model, as they must — and the weighted run used 13% fewer
likelihood calls. `mass_eep` comes out about 2% further and 0.01–0.02 dex older
in both of two seeds; reweighting its samples by the exact prior ratio (∝ dlogAge/dEEP,
which minimint returns) lands on the uniform-age answer, so that shift is the
EEP prior, not sampling.

**K giants** (MagE 18, X-shooter 33; `PARALLAX_MODE=likelihood
DISTANCE_PRIOR=flat`, 4–13 Gyr, nlive 5000, one seed; the analysis evaluates
the model at each posterior sample):

| | MagE med \|Δlog g\| | >3σ | Σ lnZ | calls | X-shooter med \|Δlog g\| | >3σ | Σ lnZ | calls |
|---|---|---|---|---|---|---|---|---|
| `mass_age` | 1.45σ | 2 | −385.0 | 48.2 M | 1.19σ | 5 | −831.0 | 91.4 M |
| `mass_eep` | 1.13σ | 2 | −300.1 | 39.3 M | 1.08σ | 3 | −723.5 | 81.2 M |
| `mass_age_weighted` | 1.22σ | **1** | −374.7 | **38.1 M** | **0.88σ** | **2** | −793.6 | **76.0 M** |

Two different effects are mixed in that table, and they separate cleanly.

*Sampling (`mass_age_weighted` against `mass_age` — same prior, same
likelihood).* The weighted runs find 10.3 (MagE) and 37.4 (X-shooter) more nats
of evidence, by more than 1 nat on 4/18 and 9/33 stars, up to 12.5;
`mass_age` is ahead by more than 1 nat on 0/18 and 3/33, by at most 2.2.
Since the model is the same, the extra evidence is posterior mass `mass_age`
did not find, and it sits at the luminous end of the track: `mass_age`
distances are a median 9% (MagE) and 2% (X-shooter) shorter, and up to 38%
shorter for individual stars. Note that lnZ is a blunt instrument here —
missing 60% of the posterior costs under 1 nat — so several MagE stars differ
by 10–25% in distance at ΔlnZ ≈ 0.

*Prior (`mass_eep` against `mass_age_weighted`).* Uniform-in-EEP does two
things to the prior, not one.

1. Inside its support it reweights phases by dEEP/dlogAge — by how fast the
   star is changing, not by how long it stays there. Relative to uniform log
   age, it weights the RGB at log g = 2.5 about 60× and at log g = 1.0 about
   550× the main sequence, a factor 2.6 across log g = 2.2 ± 0.3. Where this is
   the only difference, it is small on these data: on the 20 X-shooter stars
   whose uniform-age posterior lies inside the EEP support, the distances agree
   to a median 1.4% and none differs by 10%, and reweighting the `mass_eep`
   samples by the prior ratio brings that to 0.1%.
2. **Its support ends at `EEP_MAX = 808`, the start of the TP-AGB.** On 3/18
   MagE and 6/33 X-shooter stars, 64–100% of the uniform-age posterior lies on
   the TP-AGB or post-AGB (EEP 808–1424). `mass_eep` cannot represent those
   solutions and puts these stars at the RGB tip or early AGB instead, 11–36%
   nearer.
   Raising `EEP_MAX` does not fix this: the TP-AGB spans EEP 808–1409, 601 EEPs
   for about a million years of evolution, so a uniform-EEP prior would give it
   more weight than the whole main sequence and RGB together (EEP 202–605).

The summed lnZ gap between `mass_eep` and the uniform-age runs, +74.6 (MagE)
and +70.1 (X-shooter) nats, about 4 per star, is the Bayes factor between the
two priors — essentially the extra prior weight uniform-in-EEP gives the giant
branch over the time spent there. It is not a measure of sampling quality.
**This corrects how the 2026-09b result should be read:** of the "~85 nats in
favour of `mass_eep`" on the MagE giants, 10.3 were `mass_age` failing to
sample, and 74.6 are the change of prior.

The TP-AGB and post-AGB solutions are short-lived — the uniform-age prior
already charges them for it — and are still preferred, for stars whose
spectroscopic log g lies between −0.5 and +0.5. Whether to allow them is a
science decision, and should be made explicitly rather than inherited from
`EEP_MAX`.

Rerunning the 2026-09-15 MagE `mass_age` and `mass_eep` arrays unchanged, with
the same seed, reproduced every lnZ and every likelihood-call count exactly
(36 star fits).

**One cost.** At log age ≲ 6.6 the isochrone reaches past `M_MAX = 100`, the
proposal is truncated there and becomes close to uniform on [0.1, 100], so a
young low-mass solution occupies a thinner sliver of the unit cube than under
Salpeter. Irrelevant for old populations; worth knowing for young ones.

The output table is unchanged: column 1 is a logAge throughout, in the
checkpoint too, so nothing needs converting.

---

## 2026-09b — the EEP sampling coordinate, and a warning that was never heard

### NEW — `SAMPLE_COORD`: sample in EEP instead of age

`SAMPLE_COORD=mass_eep` samples (mass, EEP) and derives logAge from them, the
way MINESweeper does (Cargile et al. 2020). The default stays `mass_age`.

The two parameterisations describe the same model but differ in how much of the
giant branch the sampler can reach. At a fixed old age the whole RGB sits in a
narrow sliver of mass just below `getMaxMass(logAge, feh)`, and under
`mass_age` the likelihood cuts it off there. On a 0.9 Msun, [Fe/H] = −1.2
track the RGB spans 47.7% of a U(200, 808) EEP prior, against 1.8% of a
4–13 Gyr age prior. Drawing 4000 prior points and keeping those on the grid,
**none reached logg < 2 under `mass_age`, against 12.8% under `mass_eep`**,
and the lowest reachable logg was 3.09 against 0.56.

Measured on 18 distant K giants with spectroscopic logg, comparing the model
prediction at each posterior sample against the measurement:

| coordinate | median \|Δlog g\| | stars > 3σ | Σ lnZ |
|---|---|---|---|
| `mass_age` | 1.47σ | 2/18 | −385 |
| `mass_eep` | **1.15σ** | 2/18 | **−300** |

at otherwise identical priors and likelihood. A real but moderate improvement,
plus ~85 nats of evidence. *(2026-09c: most of those 85 nats are a Bayes factor
between the two priors, not better sampling — see the measurements there.)* `mass_age` is not broken — it is a worse coordinate
for evolved stars. **Prefer `mass_eep` for giants**; for dwarfs and subgiants
there is little to choose between them.

The reported posterior is unchanged: mass, logAge, [Fe/H], distance, E(B−V).
The EEP column is converted back to logAge before anything is summarised, so
percentiles, multimodality flags, corner plots and the output table see the
columns they always saw. Note that a dynesty **checkpoint** is written during
sampling and therefore still holds the raw EEP in column 1; anything reading a
checkpoint directly must convert for itself. The two are unambiguous by range
(logAge ≤ 10.2, EEP ≥ 200).

The implied age prior is whatever uniform-in-EEP induces; it is **not** uniform
in age. That is MINESweeper's choice too, and it is the honest description of
the trade: an age prior and an evolutionary-stage prior cannot both be uniform.

`EEP_MIN` / `EEP_MAX` default to 200 and 808.

### FIXED — mistfit's own warnings were silenced by its own filter

`core.py` raised the band-exclusion warning added in 2026-09 as a
`RuntimeWarning`, then installed `filterwarnings("ignore", RuntimeWarning)` at
import to suppress numpy/scipy sampling chatter. The filter swallowed the
warning, so the exclusion documented as happening "once and loudly" was in fact
silent in ordinary use.

The test covering it passed anyway, because it called
`warnings.simplefilter("always")` inside its `catch_warnings` block and so
overrode the very filter that caused the bug. There is now a second test that
leaves the module's own filters in place.

Fixed by giving mistfit its own `MistfitWarning(UserWarning)` category, exempt
from the RuntimeWarning filter. If you filter mistfit warnings by category,
catch `mistfit.core.MistfitWarning`.

### FIXED — off-grid EEP lookups returned 0.0, which is finite

`minimint.TheoryInterpolator.getLogAgeFromEEP` signals "off the grid" by
returning **0.0**, not NaN — an EEP past the end of the track, an [Fe/H]
outside the grid and a mass below it all come back as exactly 0.0. That is a
trap: 0.0 passes an `np.isfinite` guard and would be read as an age of one
year. Since MIST never produces a logAge below `LOGAGE_MIN = 5`, any
non-positive return is unambiguously the sentinel, and `logage_from_eep` now
translates it to NaN so that every caller's finiteness check means what it
says.

---

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
