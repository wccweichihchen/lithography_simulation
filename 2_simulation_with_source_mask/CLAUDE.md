# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the simulation

```bash
python sim.py --source <source.csv> --mask <mask.oas> [--n_svd 30] [--output PREFIX]
```

Example using inputs from `1_basic`:

```bash
python sim.py \
  --source ../1_basic/source_distribution.csv \
  --mask   ../1_basic/mask_lines.oas \
  --output lines
```

Outputs two PNG files: `<PREFIX>_aerial_image.png` and `<PREFIX>_linecut.png`.

Requires: `numpy`, `matplotlib`, `gdstk` (`pip install gdstk`).

## Architecture

`sim.py` is a single-file pipeline. Unlike `1_basic`, source and mask are **read from files** rather than generated in code.

### Pipeline stages

1. **`read_source_csv(filepath)`** — parses the CSV written by `1_basic/aerial_image_sim.py`. Comment lines (`#`) carry the full optical system metadata (wavelength, NA, sigma, grid N, dx). Data rows supply `fx_per_nm`, `fy_per_nm`, `intensity` for each non-zero source point, which are placed back onto the NxN frequency grid.

2. **`read_mask_oas(filepath, sys)`** — reads OASIS via `gdstk`. Each polygon's bounding box is rasterized onto the NxN grid using `row = round(y/dx)`, `col = round(x/dx)` — the inverse of the row-merge compression used in `write_oasis()`.

3. **`make_grids` / `make_pupil`** — unchanged from `1_basic`. Pupil is built with a fixed Zernike coma term `(n=3, m=1): 0.02 waves`.

4. **`abbe_simulation` / `hopkins_simulation`** — identical implementations to `1_basic`. Hopkins SVD truncation is controlled by `--n_svd` (default 30).

5. **Plotting** — `plot_results` (6-panel summary) and `plot_linecuts` (centre-row Abbe vs Hopkins overlay), both saving to files prefixed by `--output`.

### Unit convention

All lengths in **nm**, all spatial frequencies in **1/nm** — consistent with `1_basic`.

### Relationship to 1_basic

`1_basic/aerial_image_sim.py` generates masks and source in-code, then exports them. This module consumes those exports. The `OpticalSystem` dataclass and all simulation/plotting functions are kept in sync with `1_basic` by design.
