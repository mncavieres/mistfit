from astropy.table import Table
from minimint_fit import fit_stars_with_minimint

if __name__ == "__main__":
    tab = Table.read("/Users/mncavieres/Documents/2024-1/Delve/Data/xshooter_with_phot_gaia_minimint.fits")
    out = fit_stars_with_minimint(
        tab, output_path="/Users/mncavieres/Documents/2024-1/Delve/Data/xshooter_nested",
        nlive=5000, processes=8, debug=True
    )

# this is a test

