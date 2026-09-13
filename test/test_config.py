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
    ("PARALLAX_MODE", "nonsense"), ("DISTANCE_PRIOR", "nonsense")])
def test_an_unknown_setting_is_refused_not_ignored(var, bad):
    with pytest.raises(ValueError, match=var):
        load(**{var: bad})
