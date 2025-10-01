# minimint_fit.py
#!/usr/bin/env python
"""
Nested-sampling MIST + extinction fit using dynesty + minimint.

- Input: astropy Table with source_id, ≥3 photometric bands (minimint names) + their errors.
         Optional: Teff (+err), logg (+err), FEH or FEH_CAL (+err), PARALLAX (+err, +PARALLAX_ZPC), EBV (+err).
- Output: same Table with posterior summaries and multimodality flags appended.
- Debug: if debug=True, per-star folder with run/trace/corner and isochrone plots.

Author: Manuel Cavieres
Date: 2025-09-18
Requires: minimint, dynesty, astropy, numpy, scipy, matplotlib, tqdm, corner
"""

from __future__ import annotations
import os
import json
import warnings
from multiprocessing import Pool

import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use("Agg")  # headless-safe
import matplotlib.pyplot as plt

import dynesty
from dynesty import plotting as dyplot
from dynesty import DynamicNestedSampler

from scipy.stats import truncnorm, gaussian_kde
from scipy.signal import find_peaks
from astropy.table import Table

import minimint
from tqdm import tqdm

warnings.filterwarnings("ignore", category=RuntimeWarning)


# Physical / prior bounds
M_MIN, M_MAX = 0.1, 100.0                 # Msun
LOGAGE_MIN, LOGAGE_MAX = 5.0, 10.113943352306837  # yr, no need to go beyond the hubble time
FEH_MIN, FEH_MAX = -2.9, 0.9 # limits of MIST
DIST_MIN, DIST_MAX = 10.0, 2.0e5         # pc
EBV_MIN_DEFAULT, EBV_MAX_DEFAULT = 0.0, 1.5 # you might need to extend this for very reddened stars
ALPHA_IMF = 2.35 # Salpeter slope

# Supported (minimint) bands list (informative; we auto-detect from table)
MINIMINT_BANDS = set([
    'Bessell_I','PS_z','2MASS_Ks','SkyMapper_u','Gaia_RP_MAW','Gaia_RP_DR2Rev','Gaia_BP_EDR3',
    'PS_w','PS_r','SDSS_g','SkyMapper_v','Bessell_V','PS_y','SDSS_r','Bessell_U','Gaia_BP_DR2Rev',
    'WISE_W1','PS_i','Tycho_V','SDSS_i','WISE_W4','Gaia_G_EDR3','Gaia_RP_EDR3','SkyMapper_r',
    'Gaia_BP_MAWb','SDSS_z','Tycho_B','Bessell_B','DECam_i','Gaia_G_DR2Rev','DECam_Y','2MASS_J',
    'Kepler_Kp','DECam_u','GALEX_FUV','GALEX_NUV','PS_open','SDSS_u','DECam_r','SkyMapper_g',
    'SkyMapper_z','PS_g','Kepler_D51','DECam_g','Bessell_R','DECam_z','Hipparcos_Hp','WISE_W2',
    'TESS','2MASS_H','WISE_W3','SkyMapper_i','Gaia_G_MAW','Gaia_BP_MAWf'
])

# Extinction coefficients A_lambda/E(B-V)
# (Gaia triplet is handled with a color-dependent law if all three bands are present)
EXT_COEFF = {
    'Gaia_G_EDR3': 0.83627*3.1,
    'Gaia_BP_EDR3':1.08337*3.1,
    'Gaia_RP_EDR3':0.63439*3.1,
    'DECam_g':3.451,'DECam_r':2.646,'DECam_i':2.103,'DECam_z':1.575,'DECam_Y':1.515,
    'SDSS_u':4.871,'SDSS_g':3.560,'SDSS_r':2.681,'SDSS_i':2.400,'SDSS_z':1.899,
    'WISE_W1':0.0,'WISE_W2':0.0,'WISE_W3':0.0,'WISE_W4':0.0,
    'SkyMapper_u':5.017,'SkyMapper_v':4.750,'SkyMapper_g':3.651,'SkyMapper_r':3.032,
    'SkyMapper_i':2.325,'SkyMapper_z':1.790,
    '2MASS_J':0.987,'2MASS_H':0.531,'2MASS_Ks':0.164,
    # Common Johnson–Cousins/Bessell approximations (Rv=3.1-style)
    'Bessell_U':5.47,'Bessell_B':4.32,'Bessell_V':3.31,'Bessell_R':2.68,'Bessell_I':1.85,
}

# Error column suffixes to search for
ERR_SUFFIXES = ["_ERR", "_ERRMAG", "_SIG", "_E", "_error"]

# Per-process interpolator cache (keyed by sorted tuple of bands)
INTERP_CACHE = {}

# helper functions
def _get_interpolator_for_bands(bands):
    """Create/reuse a minimint.Interpolator for this exact band set (per process)."""
    key = tuple(sorted(bands))
    interp = INTERP_CACHE.get(key)
    if interp is None:
        INTERP_CACHE[key] = interp = minimint.Interpolator(list(key))
    return interp

def get_gaia_extinction(G, bp, rp, ebv, maxnit=10):
    """Color-dependent Gaia extinction (Babusieux-like)."""
    c1, c2, c3, c4, c5, c6, c7 = [0.9761, -0.1704, 0.0086, 0.0011, -0.0438, 0.0013, 0.0099]
    d1, d2, d3, d4, d5, d6, d7 = [1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043]
    e1, e2, e3, e4, e5, e6, e7 = [0.6104, -0.0170, -0.0026, -0.0017, -0.0078, 0.00005, 0.0006]
    A0 = 3.1 * ebv
    def F1(bprp): return np.poly1d([c1, c2, c3, c4][::-1])(bprp)+c5*A0+c6*A0**2+c7*bprp*A0
    def F2(bprp): return np.poly1d([d1, d2, d3, d4][::-1])(bprp)+d5*A0+d6*A0**2+d7*bprp*A0
    def F3(bprp): return np.poly1d([e1, e2, e3, e4][::-1])(bprp)+e5*A0+e6*A0**2+e7*bprp*A0
    curbp = bp - rp
    for _ in range(maxnit):
        AG  = F1(curbp) * A0
        Abp = F2(curbp) * A0
        Arp = F3(curbp) * A0
        curbp = bp - rp - Abp + Arp
    AG  = F1(curbp) * A0
    Abp = F2(curbp) * A0
    Arp = F3(curbp) * A0
    return AG, Abp, Arp

def _pick_err_col(name, cols):
    for suf in ERR_SUFFIXES:
        cand = name + suf
        if cand in cols:
            return cand
    return None

def _is_valid_number(x):
    if ma.isMaskedArray(x):
        return (x is not ma.masked) and np.isfinite(x.filled(np.nan))
    try:
        return np.isfinite(float(x))
    except Exception:
        return False

def _bands_for_row(row):
    cols = row.colnames if hasattr(row, 'colnames') else row.keys()
    usable = []
    for b in MINIMINT_BANDS:
        if b in cols:
            errc = _pick_err_col(b, cols)
            if errc and _is_valid_number(row[b]) and _is_valid_number(row[errc]) and float(row[errc]) > 0:
                usable.append(b)
    return usable

def _build_observed_from_row(row, bands):
    cols = row.colnames if hasattr(row, 'colnames') else row.keys()
    obs = {}
    # Teff / logg
    if 'Teff' in cols and _pick_err_col('Teff', cols):
        T, Terr = row['Teff'], row[_pick_err_col('Teff', cols)] + 100 # minimum 100K error for systematic floor
        if _is_valid_number(T) and _is_valid_number(Terr) and float(Terr) > 0:
            obs['Teff'] = (float(T), float(Terr))
    if 'logg' in cols and _pick_err_col('logg', cols):
        g, gerr = row['logg'], row[_pick_err_col('logg', cols)] + 0.1 # minimum 0.1 dex error for systematic floor
        if _is_valid_number(g) and _is_valid_number(gerr) and float(gerr) > 0:
            obs['logg'] = (float(g), float(gerr))
    # [Fe/H] prior
    feh_key = 'FEH_CAL' if 'FEH_CAL' in cols else ('FEH' if 'FEH' in cols else None)
    if feh_key and _pick_err_col(feh_key, cols):
        f, ferr = row[feh_key], row[_pick_err_col(feh_key, cols)] + 0.1 # minimum 0.1 dex error for systematic floor
        if _is_valid_number(f) and _is_valid_number(ferr) and float(ferr) > 0:
            obs['feh'] = (float(f), float(ferr))
    # Parallax prior (+ ZP)
    if 'PARALLAX' in cols and _pick_err_col('PARALLAX', cols):
        plx = row['PARALLAX']; perr = row[_pick_err_col('PARALLAX', cols)]
        if _is_valid_number(plx) and _is_valid_number(perr) and float(perr) > 0:
            plx = float(plx)
            if 'PARALLAX_ZPC' in cols and _is_valid_number(row['PARALLAX_ZPC']):
                plx += float(row['PARALLAX_ZPC'])
            if plx > 0:
                obs['parallax'] = (plx, float(perr))
    # EBV prior
    if 'EBV' in cols:
        ebv = row['EBV']; eerrc = _pick_err_col('EBV', cols)
        if _is_valid_number(ebv):
            if eerrc and _is_valid_number(row[eerrc]) and float(row[eerrc]) > 0:
                obs['ebv'] = (float(ebv), float(row[eerrc]))
            else:
                obs['ebv'] = (float(ebv), 0.1)
    # Photometry
    for b in bands:
        errc = _pick_err_col(b, cols)
        if errc is None:
            continue
        val, err = row[b], row[errc] + 0.1 # minimum 0.1 mag error for systematic floor
        if _is_valid_number(val) and _is_valid_number(err) and float(err) > 0:
            obs[b] = (float(val), float(err))
    return obs

def _choose_mode(obs):
    #print(obs)
    return 'spec' if ('Teff' in obs and 'logg' in obs) else 'phot'

def _summarize_samples(samples):
    p16 = np.percentile(samples, 16, axis=0)
    p50 = np.percentile(samples, 50, axis=0)
    p84 = np.percentile(samples, 84, axis=0)
    return p16, p50, p84

def _multimodal_flag(samples_1d):
    x = samples_1d
    x = x[np.isfinite(x)]
    if x.size < 200:
        return 0
    kde = gaussian_kde(x)
    xs = np.linspace(np.min(x), np.max(x), 400)
    pdf = kde(xs)
    prominence = 0.05 * np.max(pdf)
    peaks, _ = find_peaks(pdf, prominence=prominence)
    return 1 if len(peaks) > 1 else 0

def _ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path

# top level functions pickable for multiprocessing
def ptform_u5(u, obs, ebv_range,
              M_MIN_=M_MIN, M_MAX_=M_MAX,
              LOGAGE_MIN_=LOGAGE_MIN, LOGAGE_MAX_=LOGAGE_MAX,
              FEH_MIN_=FEH_MIN, FEH_MAX_=FEH_MAX,
              DIST_MIN_=DIST_MIN, DIST_MAX_=DIST_MAX,
              ALPHA_IMF_=ALPHA_IMF):
    """u~U(0,1)^5 -> (M, logAge, [Fe/H], distance_pc, E(B-V)) using available obs priors."""
    # Mass ~ Salpeter
    exp = 1.0 - ALPHA_IMF_
    m = (u[0] * (M_MAX_**exp - M_MIN_**exp) + M_MIN_**exp) ** (1.0 / exp)
    # logAge ~ uniform
    la = LOGAGE_MIN_ + u[1] * (LOGAGE_MAX_ - LOGAGE_MIN_)
    # [Fe/H]
    feh_prior = obs.get('feh')
    if feh_prior is not None:
        mu_f, sig_f = feh_prior
        a_f = (FEH_MIN_ - mu_f) / sig_f
        b_f = (FEH_MAX_ - mu_f) / sig_f
        feh = truncnorm.ppf(np.clip(u[2], 1e-12, 1-1e-12), a_f, b_f, loc=mu_f, scale=sig_f)
    else:
        feh = FEH_MIN_ + u[2] * (FEH_MAX_ - FEH_MIN_)
    # Distance
    plx_prior = obs.get('parallax')
    if plx_prior is not None:
        mu_p, sig_p = plx_prior
        a_p = (0.0 - mu_p) / sig_p
        b_p = (np.inf - mu_p) / sig_p
        plx = truncnorm.ppf(np.clip(u[3], 1e-12, 1-1e-12), a_p, b_p, loc=mu_p, scale=sig_p)
        d = 1.0e3 / plx
        d = np.clip(d, DIST_MIN_, DIST_MAX_)
    else:
        logd = np.log10(DIST_MIN_) + u[3] * (np.log10(DIST_MAX_) - np.log10(DIST_MIN_))
        d = 10**logd
    # E(B-V)
    
    EBV_MIN, EBV_MAX = ebv_range
    EBV_MAX = obs['ebv'] * 5 if ('ebv' in obs and obs['ebv'][1] > 0) else EBV_MAX_DEFAULT
    ebv = EBV_MIN + u[4] * (EBV_MAX - EBV_MIN)
    return np.array([m, la, feh, d, ebv], dtype=float)

def loglike_phot_theta(theta, obs, bands):
    """Photometry-only log-likelihood."""
    m, la, f, d, ebv = theta
    interp = _get_interpolator_for_bands(bands)
    try:
        mmax = interp.getMaxMass(la, f)
    except Exception:
        return -np.inf
    if not (M_MIN < m < mmax):
        return -np.inf

    model = interp(m, la, f)
    dm = 5*np.log10(d) - 5
    ll = 0.0
    gaia_trip = {'Gaia_G_EDR3','Gaia_BP_EDR3','Gaia_RP_EDR3'}
    use_gaia_color_ext = gaia_trip.issubset(set(obs.keys()))

    # non-Gaia or incomplete Gaia
    for band in bands:
        if band in gaia_trip and use_gaia_color_ext:
            continue
        if band not in obs:
            continue
        k = EXT_COEFF.get(band)
        if k is None:
            continue
        o, e = obs[band]
        mod = model[band]
        if not np.isfinite(mod):
            return -np.inf
        A = k * ebv
        ll += -0.5 * ((o - (mod + dm + A)) / e)**2

    # Gaia color-dependent extinction if triplet present
    if use_gaia_color_ext:
        AG, Abp, Arp = get_gaia_extinction(
            model['Gaia_G_EDR3'], model['Gaia_BP_EDR3'], model['Gaia_RP_EDR3'], ebv
        )
        for band, A in [('Gaia_G_EDR3', AG), ('Gaia_BP_EDR3', Abp), ('Gaia_RP_EDR3', Arp)]:
            o, e = obs[band]
            ll += -0.5 * ((o - (model[band] + dm + A)) / e)**2

    return float(ll)

def loglike_spec_theta(theta, obs, bands):
    """Photometry + spectroscopy log-likelihood."""
    ll = loglike_phot_theta(theta, obs, bands)
    if not np.isfinite(ll):
        return -np.inf
    m, la, f, d, ebv = theta
    interp = _get_interpolator_for_bands(bands)
    model = interp(m, la, f)
    if 'Teff' in obs:
        oT, eT = obs['Teff']
        ll += -0.5 * ((oT - 10**model['logteff']) / eT)**2
    if 'logg' in obs:
        og, eg = obs['logg']
        ll += -0.5 * ((og - model['logg']) / eg)**2
    return float(ll)

# main function
def fit_stars_with_minimint(
    table: Table,
    output_path: str,
    nlive: int = 3000,
    dlogz: float = 0.01,
    processes: int = 4,
    debug: bool = False,
    random_seed: int = 42,
    ebv_range: tuple[float, float] = (EBV_MIN_DEFAULT, EBV_MAX_DEFAULT),
) -> Table:
    """
    Run dynesty (internal pool) per star and append posterior summaries to the input table.

    Columns appended:
      mass_p16/p50/p84, mass_multimodal
      logage_p16/p50/p84, logage_multimodal
      feh_p16/p50/p84, feh_multimodal
      dist_pc_p16/p50/p84, dist_multimodal
      ebv_p16/p50/p84, ebv_multimodal
      lnZ, lnZ_err
    """
    rng = np.random.default_rng(random_seed)
    _ensure_dir(output_path)

    # Prepare output columns if missing
    outcols = [
        'mass_p16','mass_p50','mass_p84','mass_multimodal',
        'logage_p16','logage_p50','logage_p84','logage_multimodal',
        'feh_p16','feh_p50','feh_p84','feh_multimodal',
        'dist_pc_p16','dist_pc_p50','dist_pc_p84','dist_multimodal',
        'ebv_p16','ebv_p50','ebv_p84','ebv_multimodal',
        'lnZ','lnZ_err', 'fit_mode'
    ]
    for c in outcols:
        if c not in table.colnames:
            table[c] = np.full(len(table), np.nan)
    # make the fit_mode column string type
    if 'fit_mode' in table.colnames:
        table['fit_mode'] = table['fit_mode'].astype(str)

    for i, row in tqdm(enumerate(table), total=len(table), desc="Fitting stars"):
        sid = str(row['source_id']) if 'source_id' in row.colnames else f"row{i}"

        # Detect usable bands
        bands = _bands_for_row(row)
        if len(bands) < 3:
            # Not enough data; leave NaNs
            continue

        obs = _build_observed_from_row(row, bands)

        # EBV range per star (if prior provided)
        if 'ebv' in obs and obs['ebv'][1] > 0:
            ebv_min = max(0.0, obs['ebv'][0] - 5*obs['ebv'][1])
            ebv_max = min(EBV_MAX_DEFAULT, obs['ebv'][0] + 5*obs['ebv'][1])
            ebv_range_row = (ebv_min, max(ebv_min + 1e-3, ebv_max))
        else:
            ebv_range_row = ebv_range

        # Mode & loglike
        mode = _choose_mode(obs)
        print(f"[{i+1}/{len(table)}] Fitting {sid} with {mode}-only mode, {len(bands)} bands")
        loglike = loglike_spec_theta if mode == 'spec' else loglike_phot_theta
        table['fit_mode'][i] = mode

        # Args (must be picklable)
        logl_args = (obs, bands)
        ptform_args = (obs, ebv_range_row)

        ndim = 5
        # Use dynesty's internal multiprocessing pool (picklable callables!)
        with Pool(processes=processes) as pool_in:
            sampler = DynamicNestedSampler(
                loglike, ptform_u5, ndim,
                bound='multi', sample='rwalk',
                pool=pool_in, queue_size=processes,
                logl_args=logl_args, ptform_args=ptform_args,
                rstate=rng,
            )
            ckpt_file = os.path.join(output_path, f'{sid}_checkpoint.h5') if debug else None
            sampler.run_nested(
                nlive_init=nlive,
                dlogz_init=dlogz,
                #print_progress=debug,
                checkpoint_file=ckpt_file
            )
        res = sampler.results

        # Equal-weight posterior samples
        samples = res.samples_equal()  # (Nsamp, 5)
        if samples.size == 0:
            continue

        # Summaries
        p16, p50, p84 = _summarize_samples(samples)
        (m16, la16, f16, d16, e16) = p16
        (m50, la50, f50, d50, e50) = p50
        (m84, la84, f84, d84, e84) = p84

        # Multimodality flags
        mflag  = _multimodal_flag(samples[:,0])
        laflag = _multimodal_flag(samples[:,1])
        fflag  = _multimodal_flag(samples[:,2])
        dflag  = _multimodal_flag(samples[:,3])
        eflag  = _multimodal_flag(samples[:,4])

        # Write back
        table['mass_p16'][i] = m16; table['mass_p50'][i] = m50; table['mass_p84'][i] = m84
        table['mass_multimodal'][i] = mflag

        table['logage_p16'][i] = la16; table['logage_p50'][i] = la50; table['logage_p84'][i] = la84
        table['logage_multimodal'][i] = laflag

        table['feh_p16'][i] = f16; table['feh_p50'][i] = f50; table['feh_p84'][i] = f84
        table['feh_multimodal'][i] = fflag

        table['dist_pc_p16'][i] = d16; table['dist_pc_p50'][i] = d50; table['dist_pc_p84'][i] = d84
        table['dist_multimodal'][i] = dflag

        table['ebv_p16'][i] = e16; table['ebv_p50'][i] = e50; table['ebv_p84'][i] = e84
        table['ebv_multimodal'][i] = eflag

        table['lnZ'][i] = res.logz[-1] if len(res.logz) else np.nan
        table['lnZ_err'][i] = res.logzerr[-1] if len(res.logzerr) else np.nan

        # Debug plots
        if debug:
            star_dir = _ensure_dir(os.path.join(output_path, sid))
            # Run plots
            try:
                rfig, _ = dyplot.runplot(res)
                rfig.savefig(os.path.join(star_dir, 'runplot.png'), dpi=150, bbox_inches='tight')
                plt.close(rfig)
            except Exception:
                pass
            # Trace plots
            try:
                tfig, _ = dyplot.traceplot(res)
                tfig.savefig(os.path.join(star_dir, 'traceplot.png'), dpi=150, bbox_inches='tight')
                plt.close(tfig)
            except Exception:
                pass
            # Corner plot
            try:
                import corner
                labels = ['Mass','log10(age/yr)','[Fe/H]','Distance','E(B-V)']
                cfig = corner.corner(samples, labels=labels, bins=30, smooth=1.5,
                                     show_titles=True, title_kwargs={'fontsize':9})
                cfig.savefig(os.path.join(star_dir, 'corner.png'), dpi=150, bbox_inches='tight')
                plt.close(cfig)
            except Exception:
                pass
            # Isochrone check (requires Gaia G + Teff)
            try:
                has_gaiaG = ('Gaia_G_EDR3' in bands) and ('Gaia_G_EDR3' in obs)
                has_teff  = 'Teff' in obs
                if has_gaiaG and has_teff:
                    best_mass, best_logage, best_feh, best_dist, best_ebv = p50
                    interp = _get_interpolator_for_bands(bands)
                    masses_grid = np.linspace(M_MIN, interp.getMaxMass(best_logage, best_feh), 200)
                    iso = interp(masses_grid, best_logage, best_feh)
                    teffs = 10**iso['logteff']
                    dm = 5*np.log10(best_dist) - 5
                    kG = EXT_COEFF['Gaia_G_EDR3']
                    Gs = iso['Gaia_G_EDR3'] + dm + kG * best_ebv
                    obs_T, err_T = obs['Teff']
                    obs_G, err_G = obs['Gaia_G_EDR3']
                    fig, ax = plt.subplots(figsize=(6,8))
                    ax.plot(teffs, Gs, '-', lw=2, label='Best-fit isochrone')
                    ax.errorbar(obs_T, obs_G, xerr=err_T, yerr=err_G,
                                fmt='o', ms=5, color='black', label='Observed')
                    ax.invert_xaxis(); ax.invert_yaxis()
                    ax.set_xlabel(r'$T_{\rm eff}$ [K]')
                    ax.set_ylabel('Gaia G [mag]')
                    ax.set_title('Isochrone check')
                    ax.set_xlim(max(3000, teffs.min()), min(10000, teffs.max()))
                    ax.legend()
                    fig.savefig(os.path.join(star_dir, 'isochrone_check.png'), dpi=150, bbox_inches='tight')
                    plt.close(fig)
            except Exception:
                pass
            # JSON summary
            try:
                with open(os.path.join(star_dir, 'summary.json'), 'w') as f:
                    json.dump({
                        'mode': mode,
                        'bands': bands,
                        'p50': {'mass': float(m50), 'logage': float(la50), 'feh': float(f50),
                                'dist_pc': float(d50), 'ebv': float(e50)},
                        'lnZ': float(table['lnZ'][i]),
                        'lnZ_err': float(table['lnZ_err'][i])
                    }, f, indent=2)
            except Exception:
                pass

        # Save updated table
        out_file = os.path.join(output_path, 'fit_results.fits')
        try:
            table.write(out_file, overwrite=True)
        except Exception:
            # fall back to ECSV if FITS column types collide
            table.write(out_file.replace('.fits', '.ecsv'), overwrite=True)
    return table


# if __name__ == "__main__":
#     import argparse
#     ap = argparse.ArgumentParser(description="Run minimint nested fits on an input catalog.")
#     ap.add_argument("input_table", help="Path to input FITS/ECSV/etc.")
#     ap.add_argument("--outdir", required=True, help="Output directory")
#     ap.add_argument("--nlive", type=int, default=3000)
#     ap.add_argument("--dlogz", type=float, default=0.01)
#     ap.add_argument("--procs", type=int, default=4)
#     ap.add_argument("--debug", action="store_true")
#     args = ap.parse_args()

#     tab = Table.read(args.input_table)
#     _ = fit_stars_with_minimint(
#         tab, output_path=args.outdir,
#         nlive=args.nlive, dlogz=args.dlogz,
#         processes=args.procs, debug=args.debug
#     )

# if __name__ == "__main__":
#     from astropy.table import Table

#     data = Table.read('/Users/mncavieres/Documents/2024-1/Delve/Data/xshooter_with_phot_gaia_minimint.fits')
#     out = fit_stars_with_minimint(
#     data,
#     output_path="/Users/mncavieres/Documents/2024-1/Delve/Data/xshooter_nested",
#     nlive=5000,
#     processes=8,
#     debug=True
# )
# print(out.colnames[-20:])  # see the appended posterior columns