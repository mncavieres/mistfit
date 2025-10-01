# setup.py
from pathlib import Path
from setuptools import setup, find_packages

BASE_DIR = Path(__file__).parent
readme_path = BASE_DIR / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""

setup(
    name="mistfit",
    version="0.1.0",  # switch to setuptools_scm later if you want tag-based versions
    description="Fit MIST models using MiniMINT with dynesty nested sampling (photometry + spectroscopic constraints).",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Manuel Cavieres",
    license="MIT",
    url="https://github.com/mncavieres/MIST-stellar-fitter",
    project_urls={
        "Issues": "https://github.com/mncavieres/MIST-stellar-fitter/issues",
    },
    keywords=["MIST", "stellar", "astrophysics", "dynesty", "minimint", "isochrones"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    python_requires=">=3.9",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,  # uses MANIFEST.in for sdists
    install_requires=[
        "numpy>=1.23",
        "scipy>=1.10",
        "astropy>=5",
        "dynesty>=2.1",
        "minimint>=0.5",      # adjust to the version you need
        "matplotlib>=3.6",
        "tqdm>=4.64",
        "corner>=2.2",       # keep optional if you only use it in debug
    ],
    extras_require={
        "dev": ["pytest", "pytest-cov", "ruff", "mypy", "build", "twine"],
        "plot": ["corner>=2.2"],
    },
    entry_points={
        "console_scripts": [
            "mist-fit=mistfit.cli:main",
        ]
    },
)
