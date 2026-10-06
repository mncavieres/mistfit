"""Pin the behaviour of the 2026-09 fixes.

Each test corresponds to a claim in CHANGELOG.md, so a regression shows up as a
failing assertion rather than as a quietly different number months later.

Nothing here touches the sampler or the isochrone grids: these exercise
configuration handling and pure-python helpers, and run in seconds.
"""
import importlib
import os
import sys
import warnings

import numpy as np
import pytest

CONFIG_VARS = (
    "PARALLAX_MODE", "DISTANCE_PRIOR", "RNG_PER_STAR",
    "TEFF_ERR_ADDITIVE", "LOGG_ERR_ADDITIVE", "FEH_ERR_ADDITIVE",
    "PHOT_ERR_ADDITIVE", "TEFF_ERR_FLOOR", "LOGG_ERR_FLOOR",
    "FEH_ERR_FLOOR", "PHOT_ERR_FLOOR",
    "SAMPLE_COORD", "EEP_MIN", "EEP_MAX",
)


def load(**env):
    """Re-import mistfit.core with these environment settings."""
    for k in CONFIG_VARS:
        os.environ.pop(k, None)
    os.environ.update({k: str(v) for k, v in env.items()})
    sys.modules.pop("mistfit.core", None)
    sys.modules.pop("mistfit", None)
    return importlib.import_module("mistfit.core")


@pytest.fixture(autouse=True)
def _clean_env():
    saved = {k: os.environ.get(k) for k in CONFIG_VARS}
    yield
    for k, v in saved.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v
    sys.modules.pop("mistfit.core", None)
    sys.modules.pop("mistfit", None)


class Row(dict):
    """Stand-in for an astropy Row."""
    @property
    def colnames(self):
        return list(self)


# --------------------------------------------------------------- band ordering
def test_band_list_is_sorted_not_hash_ordered():
    m = load()
    row = Row()
    for b in ("DECam_z", "DECam_g", "DECam_r", "DECam_i"):
        row[b], row[b + "_ERR"] = 15.0, 0.01
    assert m._bands_for_row(row) == sorted(m._bands_for_row(row))


def test_band_list_is_stable_across_reimports():
    assert load()._bands_for_row.__module__ == "mistfit.core"
    first = sorted(load().MINIMINT_BANDS)
    assert first == sorted(load().MINIMINT_BANDS)


# ------------------------------------------------------- extinction coverage
def test_vista_bands_present_with_their_own_coefficients():
    m = load()
    vista = sorted(b for b in m.MINIMINT_BANDS if b.startswith("VISTA_"))
    assert vista == ["VISTA_H", "VISTA_J", "VISTA_Ks", "VISTA_Y"]
    assert all(b in m.EXT_COEFF for b in vista)
    # VISTA is not 2MASS: A_Ks/E(B-V) differs by a factor 2.4.
    assert m.EXT_COEFF["VISTA_Ks"] != m.EXT_COEFF["2MASS_Ks"]


def test_band_without_extinction_coefficient_is_excluded_and_warned():
    m = load()
    orphan = sorted(b for b in m.MINIMINT_BANDS
                    if b not in m.EXT_COEFF and b not in m._GAIA_TRIPLET)
    if not orphan:
        pytest.skip("every band has a coefficient")
    row = Row()
    for b in orphan[:2] + ["DECam_g", "DECam_r", "DECam_i"]:
        row[b], row[b + "_ERR"] = 15.0, 0.01
    m._WARNED_NO_EXT_COEFF.clear()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        bands = m._bands_for_row(row)
    assert not set(bands) & set(orphan), "coefficient-less band leaked into the fit"
    assert {"DECam_g", "DECam_r", "DECam_i"} <= set(bands)
    assert any("EXT_COEFF" in str(c.message) for c in caught), "exclusion was silent"


# --------------------------------------------------------------- distance prior
@pytest.mark.parametrize("mode", ["loguniform", "flat", "volume"])
def test_distance_prior_is_a_valid_inverse_cdf(mode):
    m = load(DISTANCE_PRIOR=mode)
    lo, hi = 10.0, 2.0e5
    edges = m._draw_distance(np.array([0.0, 1.0]), lo, hi)
    assert edges[0] == pytest.approx(lo, rel=1e-12)
    assert edges[1] == pytest.approx(hi, rel=1e-12)
    u = np.linspace(0, 1, 5000)
    d = m._draw_distance(u, lo, hi)
    assert np.all(np.diff(d) >= 0)
    assert d.min() >= lo - 1e-6 and d.max() <= hi + 1e-6


def test_distance_priors_are_ordered_by_how_much_near_space_they_favour():
    lo, hi = 10.0, 2.0e5
    u = (np.arange(200001) + 0.5) / 200001
    frac = {}
    for mode in ("loguniform", "flat", "volume"):
        d = load(DISTANCE_PRIOR=mode)._draw_distance(u, lo, hi)
        frac[mode] = np.mean(d < 1.0e4)
    assert frac["loguniform"] > frac["flat"] > frac["volume"]
    # the historical default puts most of its mass inside 10 kpc
    assert frac["loguniform"] > 0.5
    assert frac["volume"] < 0.01


def test_loguniform_remains_the_default():
    assert load().DISTANCE_PRIOR == "loguniform"


# --------------------------------------------------------------- parallax mode
@pytest.mark.parametrize("mode,keeps_negative", [
    ("published", False), ("prior", True), ("likelihood", True)])
def test_negative_parallax_is_kept_unless_reproducing_old_behaviour(mode, keeps_negative):
    m = load(PARALLAX_MODE=mode)
    neg = Row(PARALLAX=-0.05, PARALLAX_ERR=0.06)
    pos = Row(PARALLAX=+0.05, PARALLAX_ERR=0.06)
    assert ("parallax" in m._build_observed_from_row(neg, [])) is keeps_negative
    assert "parallax" in m._build_observed_from_row(pos, [])


def test_prior_is_the_default_parallax_mode():
    assert load().PARALLAX_MODE == "prior"


def test_parallax_enters_prior_or_likelihood_but_never_both():
    """The prior transform must use the parallax in "prior" mode and ignore it
    in "likelihood" mode, where the measurement belongs to the likelihood."""
    from scipy.stats import truncnorm

    mu, sig = 0.02, 0.05
    obs = {"parallax": (mu, sig)}
    u = np.full(5, 0.5)

    d_prior = load(PARALLAX_MODE="prior").ptform_u5(u, obs, (0.0, 1.5))[3]
    # 1e3 / the median of TruncNorm(mu, sig, lower=0). Note that is not 1e3/mu:
    # truncating away the negative tail shifts the median up.
    expected = 1e3 / truncnorm.ppf(0.5, (0.0 - mu) / sig, np.inf, loc=mu, scale=sig)
    assert d_prior == pytest.approx(expected, rel=1e-9)

    m = load(PARALLAX_MODE="likelihood")
    d_lik = m.ptform_u5(u, obs, (0.0, 1.5))[3]
    # ignores the parallax, so it must be the DISTANCE_PRIOR median instead
    assert d_lik == pytest.approx(
        m._draw_distance(0.5, m.DIST_MIN, m.DIST_MAX), rel=1e-9)
    assert d_lik != pytest.approx(d_prior, rel=1e-6)


# ---------------------------------------------------------------- error budget
def test_additive_defaults_reproduce_the_historical_error_budget():
    m = load()
    assert (m.TEFF_ERR_ADDITIVE, m.LOGG_ERR_ADDITIVE,
            m.FEH_ERR_ADDITIVE, m.PHOT_ERR_ADDITIVE) == (100.0, 0.1, 0.1, 0.1)
    assert all(f is None for f in (m.TEFF_ERR_FLOOR, m.LOGG_ERR_FLOOR,
                                   m.FEH_ERR_FLOOR, m.PHOT_ERR_FLOOR))
    assert m._inflate(0.02, 0.1, None) == pytest.approx(0.12)


def test_a_floor_lifts_an_implausibly_small_error_and_leaves_a_large_one():
    m = load(LOGG_ERR_FLOOR=0.3)
    assert m._inflate(0.0004, 0.1, m.LOGG_ERR_FLOOR) == pytest.approx(0.3)
    assert m._inflate(0.5, 0.1, m.LOGG_ERR_FLOOR) == pytest.approx(0.6)


# ---------------------------------------------------------------------- seeding
def test_per_star_seeding_is_the_default():
    assert load().RNG_PER_STAR is True
    assert load(RNG_PER_STAR=0).RNG_PER_STAR is False


def test_a_stars_random_stream_depends_only_on_which_star_it_is():
    def draws(seed, key):
        return np.random.default_rng(
            np.random.SeedSequence([int(seed), int(key)])).random(4)
    assert np.array_equal(draws(42, 12345), draws(42, 12345))
    assert not np.array_equal(draws(42, 12345), draws(42, 12346))
    assert not np.array_equal(draws(42, 12345), draws(43, 12345))


# ------------------------------------------------------------------- validation
@pytest.mark.parametrize("var,bad", [
    ("PARALLAX_MODE", "nonsense"), ("DISTANCE_PRIOR", "nonsense"),
    ("SAMPLE_COORD", "nonsense")])
def test_an_unknown_setting_is_refused_not_ignored(var, bad):
    with pytest.raises(ValueError, match=var):
        load(**{var: bad})


# ------------------------------------------------------------ warning plumbing
def test_mistfit_warnings_survive_the_modules_own_runtime_filter():
    """The band-exclusion warning must reach the user under real conditions.

    The test above catches it with simplefilter("always"), which overrides the
    module's own filters -- so it passed even while the warning was, in
    practice, silent: core.py raised it as a RuntimeWarning and then installed
    filterwarnings("ignore", RuntimeWarning) at import. Here the module's own
    filters are left in place, which is what a user actually gets.
    """
    m = load()
    orphan = sorted(b for b in m.MINIMINT_BANDS
                    if b not in m.EXT_COEFF and b not in m._GAIA_TRIPLET)
    if not orphan:
        pytest.skip("every band has a coefficient")
    row = Row()
    for b in orphan[:1] + ["DECam_g", "DECam_r", "DECam_i"]:
        row[b], row[b + "_ERR"] = 15.0, 0.01
    m._WARNED_NO_EXT_COEFF.clear()
    with warnings.catch_warnings(record=True) as caught:
        # deliberately NO simplefilter: keep whatever core.py installed
        m._bands_for_row(row)
    assert any("EXT_COEFF" in str(c.message) for c in caught), (
        "band exclusion was silent under the module's own warning filters")


def test_numpy_runtime_noise_is_still_suppressed():
    """Silencing sampler chatter is the reason the filter exists; keep it."""
    m = load()
    with warnings.catch_warnings(record=True) as caught:
        warnings.warn("overflow encountered in exp", RuntimeWarning)
    assert not caught, "RuntimeWarning noise is no longer suppressed"
    assert issubclass(m.MistfitWarning, UserWarning)


# ------------------------------------------------------------ sampling coordinate
def test_mass_age_remains_the_default_coordinate():
    assert load().SAMPLE_COORD == "mass_age"
    assert load(SAMPLE_COORD="mass_eep").SAMPLE_COORD == "mass_eep"
    assert load(SAMPLE_COORD="mass_age_weighted").SAMPLE_COORD == "mass_age_weighted"


def test_second_coordinate_is_an_eep_only_under_mass_eep():
    """ptform's second element switches meaning; its range is the giveaway."""
    u = np.array([0.5, 0.5, 0.5, 0.5, 0.5])
    obs, ebv_range = {}, (0.0, 0.5)

    age = load(SAMPLE_COORD="mass_age").ptform_u5(u, obs, ebv_range)[1]
    assert 5.0 <= age <= 10.2, "logAge left its prior range"

    m = load(SAMPLE_COORD="mass_eep")
    eep = m.ptform_u5(u, obs, ebv_range)[1]
    assert m.EEP_MIN <= eep <= m.EEP_MAX
    assert eep > 100, "second coordinate is not an EEP under mass_eep"


def test_eep_bounds_are_configurable():
    m = load(SAMPLE_COORD="mass_eep", EEP_MIN=300, EEP_MAX=400)
    assert (m.EEP_MIN, m.EEP_MAX) == (300.0, 400.0)
    eep = m.ptform_u5(np.full(5, 0.5), {}, (0.0, 0.5))[1]
    assert 300.0 <= eep <= 400.0


@pytest.mark.parametrize("mass,eep,feh", [
    (0.9, 2000, -1.2),     # EEP past the end of the track
    (0.9, 600, -5.0),      # [Fe/H] off the grid
    (0.05, 600, -1.2),     # mass below the grid
    (1e6, 1e6, 99.0),      # nothing about this is on the grid
])
def test_off_grid_points_map_to_nan_not_to_zero(mass, eep, feh):
    """minimint returns 0.0 off-grid; we must not pass that on as an age.

    0.0 is finite, so an isfinite() guard lets it through, and it would then
    be read as an age of one year. logage_from_eep translates the sentinel so
    that every caller's finiteness check means what it says.
    """
    m = load(SAMPLE_COORD="mass_eep")
    assert np.isnan(m.logage_from_eep(mass, eep, feh))


def test_on_grid_points_still_return_a_real_age():
    m = load(SAMPLE_COORD="mass_eep")
    rgb = m.logage_from_eep(0.9, 600, -1.2)
    ms = m.logage_from_eep(0.9, 250, -1.2)
    assert m.LOGAGE_MIN <= ms < rgb <= m.LOGAGE_MAX, (ms, rgb)


# ------------------------------------------------- weighted mass proposal
@pytest.mark.parametrize("iso_maxmass,mmax", [
    (0.865, 100.0),    # old metal-poor isochrone: the case the proposal is for
    (19.8, 100.0),     # young, but still inside the mass prior
    (300.0, 100.0),    # youngest: the isochrone runs past M_MAX -> truncated
    (0.12, 100.0),     # barely above M_MIN: the break points collapse onto it
])
def test_giant_branch_proposal_is_a_normalised_inverse_cdf(iso_maxmass, mmax):
    """pdf must be the density of ppf, or the weight ln S - ln q is wrong."""
    m = load(SAMPLE_COORD="mass_age_weighted")
    q = m.GiantBranchMassProposal(iso_maxmass, 0.1, mmax)
    top = min(iso_maxmass, mmax)
    u = np.linspace(0, 1, 1001)
    x = np.array([q.ppf(v) for v in u])
    assert x[0] == pytest.approx(0.1) and x[-1] == pytest.approx(top)
    assert np.all(np.diff(x) >= 0)
    assert max(abs(q.cdf(xx) - v) for xx, v in zip(x, u)) < 1e-9
    # pdf is the derivative of cdf, away from the break points
    for xx in np.linspace(0.1, top, 23)[1:-1]:
        h = 1e-7 * top
        if min(abs(xx - b) for b in (q.m1, q.m2)) < 10 * h:
            continue
        assert q.pdf(xx) == pytest.approx((q.cdf(xx + h) - q.cdf(xx - h)) / (2 * h), rel=1e-4)
    assert q.pdf(0.1) == 0 and q.pdf(top) == 0 and q.pdf(top * 1.01) == 0


def test_giant_branch_proposal_crowds_the_isochrone_tip():
    """The point of it: far more draws near maxMass than Salpeter gives."""
    m = load(SAMPLE_COORD="mass_age_weighted")
    mx = 0.865
    q = m.GiantBranchMassProposal(mx, m.M_MIN, m.M_MAX)
    u = (np.arange(20000) + 0.5) / 20000
    frac_q = np.mean([q.ppf(v) > q.m1 for v in u])
    e = 1 - m.ALPHA_IMF
    S = lambda x: (x**e - m.M_MIN**e) / (m.M_MAX**e - m.M_MIN**e)
    frac_salpeter = S(mx) - S(q.m1)
    assert frac_q > 0.2 and frac_salpeter < 0.002
    assert frac_q / frac_salpeter > 100
    # and ~10% in the last 0.03%, the TP-AGB / post-AGB: what lets the weighted
    # runs find modes there that mass_age misses (CHANGELOG 2026-09c)
    tip_q = 1.0 - q.c2
    tip_salpeter = S(mx) - S(q.m2)
    assert q.m2 == pytest.approx(mx * (1 - 10**-3.51))
    assert tip_q > 0.09 and tip_q / tip_salpeter > 1000


def test_weight_undoes_the_proposal_exactly():
    """q(m) * exp(weight) must be the Salpeter density mass_age samples from,
    normalisation included -- that is what makes lnZ comparable to mass_age."""
    from scipy.integrate import quad

    m = load(SAMPLE_COORD="mass_age_weighted")
    m._max_mass = lambda la, feh: 0.865        # no grid needed
    q = m.GiantBranchMassProposal(0.865, m.M_MIN, m.M_MAX)
    e = 1 - m.ALPHA_IMF
    for x in (0.15, 0.5, 0.86, 0.8649):
        salpeter = e * x**-m.ALPHA_IMF / (m.M_MAX**e - m.M_MIN**e)
        assert q.pdf(x) * np.exp(m.mass_weight_log(x, 10.0, -1.2)) == pytest.approx(salpeter, rel=1e-10)
    # integrated, the weighted proposal carries exactly mass_age's on-isochrone prior volume
    brk = [q.m1, q.m2]
    got = quad(lambda x: q.pdf(x) * np.exp(m.mass_weight_log(x, 10.0, -1.2)),
               m.M_MIN, 0.865, points=brk, limit=200)[0]
    want = (0.865**e - m.M_MIN**e) / (m.M_MAX**e - m.M_MIN**e)
    assert got == pytest.approx(want, rel=1e-8)


def test_weight_rejects_what_mass_age_rejects():
    m = load(SAMPLE_COORD="mass_age_weighted")
    m._max_mass = lambda la, feh: 0.865
    assert m.mass_weight_log(0.9, 10.0, -1.2) == -np.inf       # past the tip
    assert m.mass_weight_log(0.1, 10.0, -1.2) == -np.inf       # at M_MIN
    m._max_mass = lambda la, feh: float("nan")                 # off the grid
    assert m.mass_weight_log(0.5, 10.0, -1.2) == -np.inf


def test_weighted_ptform_draws_age_as_mass_age_and_mass_from_the_proposal():
    u = np.array([0.9, 0.3, 0.5, 0.5, 0.5])
    plain = load(SAMPLE_COORD="mass_age").ptform_u5(u, {}, (0.0, 0.5))
    m = load(SAMPLE_COORD="mass_age_weighted")
    m._max_mass = lambda la, feh: 0.865
    weighted = m.ptform_u5(u, {}, (0.0, 0.5))
    # every coordinate but the mass is untouched
    assert np.allclose(weighted[1:], plain[1:])
    assert weighted[0] == pytest.approx(
        m.GiantBranchMassProposal(0.865, m.M_MIN, m.M_MAX).ppf(0.9))
    assert m.M_MIN < weighted[0] < 0.865


def test_proposal_and_weight_see_the_same_mass_bounds_after_an_override():
    """ptform_u5's defaults are frozen at import; the weighted proposal must not
    use them, or overriding core.M_MAX would change q in the likelihood only."""
    m = load(SAMPLE_COORD="mass_age_weighted")
    m._max_mass = lambda la, feh: 300.0        # a young isochrone past M_MAX
    m.M_MAX = 50.0
    u = np.array([0.999, 0.5, 0.5, 0.5, 0.5])
    mass = m.ptform_u5(u, {}, (0.0, 0.5))[0]
    assert mass <= 50.0
    assert np.isfinite(m.mass_weight_log(mass * 0.999, 6.0, 0.0))
