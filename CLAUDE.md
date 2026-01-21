# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

FRB is a Python package for Fast Radio Burst (FRB) calculations, observations, estimations, and analysis. It provides tools for DM (Dispersion Measure) calculations, host galaxy analysis, FRB-galaxy associations, and survey data handling.

## Build and Test Commands

```bash
# Install dependencies
pip install -r frb/requirements.txt
pip install git+https://github.com/FRBs/ne2001.git#egg=ne2001
pip install git+https://github.com/FRBs/astropath.git#egg=astropath
pip install git+https://github.com/linetools/linetools#egg=linetools

# Install package in development mode
pip install -e .

# Run all tests
pytest --pyargs frb

# Run specific test file
pytest frb/tests/test_frb.py

# Run with tox (specific Python version)
tox -e 3.11-test

# Run with coverage
tox -e 3.11-test-cov
```

## Architecture

### Core Classes

- **`frb.frb.FRB`**: Main class representing an observed FRB event. Load by name with `FRB.by_name('FRB20121102A')`. Contains coordinates, DM, RM, redshift, and links to host galaxies.

- **`frb.galaxies.frbgalaxy.FRBGalaxy`** / **`FRBHost`**: Classes for galaxies in FRB fields and confirmed host galaxies. Store photometry, redshifts, morphology, and derived properties.

- **`frb.associate.frbassociate.FRBAssociate`**: Extends the astropath PATH class for probabilistic FRB-galaxy association analysis.

### Key Modules

- **`frb/dm/`**: Dispersion measure calculations
  - `igm.py`: IGM DM contributions, cosmic DM estimates
  - `host.py`: Host galaxy DM contributions
  - `cosmic.py`: Cosmological DM calculations
  - `prob_dmz.py`: DM-redshift probability distributions

- **`frb/galaxies/`**: Host galaxy analysis
  - `photom.py`: Photometry handling and conversions
  - `nebular.py`: Nebular emission line analysis
  - `eazy.py` / `cigale.py`: SED fitting interfaces

- **`frb/surveys/`**: Survey data access (PanSTARRS, SDSS, DECaLS, WISE, 2MASS, etc.)

- **`frb/halos/`**: Halo mass function and foreground halo DM models

### Data Storage

FRB and host galaxy data are stored as JSON files in `frb/data/FRBs/` and `frb/data/Galaxies/`. Access via class methods like `FRB.by_name()` or `FRBHost.by_frb()`.

### Cosmology

Default cosmology is Planck18, defined in `frb/defs.py` as `frb_cosmo`. All calculations use this unless explicitly overridden.

### External Dependencies

Critical external packages (not on PyPI):
- `ne2001`: Galactic electron density model for ISM DM
- `astropath`: PATH probabilistic association framework
- `linetools`: Spectroscopic analysis utilities
