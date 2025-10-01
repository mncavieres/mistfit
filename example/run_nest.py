from astropy.table import Table
#from core import fit_stars_with_minimint
from mistfit import fit_stars_with_minimint

if __name__ == "__main__":
    tab = Table.read("photometry_spec_table.fits")
    out = fit_stars_with_minimint(
        tab, output_path="example",
        nlive=5000, processes=8, debug=True
    )

# this is a test

