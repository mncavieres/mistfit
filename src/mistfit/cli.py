# src/mist_stellar_fitter/cli.py
import argparse
from astropy.table import Table
from .core import fit_stars_with_minimint

def main():
    ap = argparse.ArgumentParser(
        prog="mist-fit",
        description="Nested-sampling MIST + extinction fit using dynesty + minimint."
    )
    ap.add_argument("input_table", help="Path to input FITS/ECSV/etc.")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--nlive", type=int, default=3000)
    ap.add_argument("--dlogz", type=float, default=0.01)
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    tab = Table.read(args.input_table)
    _ = fit_stars_with_minimint(
        tab,
        output_path=args.outdir,
        nlive=args.nlive,
        dlogz=args.dlogz,
        processes=args.procs,
        debug=args.debug,
    )

if __name__ == "__main__":
    main()
