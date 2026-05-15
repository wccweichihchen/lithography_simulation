# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Run aerial image simulation:**
```bash
python sim.py --source <source.csv> --mask <mask.oas> [--n_svd 30] [--linecut_y Y_NM] [--output PREFIX]
```

**Generate a custom layout:**
```bash
python generate_layout.py [--seed SEED] [--output PATH]
```

**End-to-end example using the bundled source and generated layout:**
```bash
python generate_layout.py --output multi_L-shaped.oas
python sim.py --source source_distribution.csv --mask multi_L-shaped.oas --output multi_L-shaped --linecut_y 100
```

Requires: `numpy`, `matplotlib`, `gdstk` (`pip install gdstk`).

## Architecture

Two scripts; source and mask are always **read from files**, never generated in code.

### `sim.py` — simulation pipeline

1. **`read_source_csv(filepath)`** — parses the CSV from `1_basic`. Comment header lines (`#`) carry optical system metadata (wavelength, NA, sigma, grid N, dx), which are used to reconstruct the `OpticalSystem`. Data rows supply `fx_per_nm`, `fy_per_nm`, `intensity` for each non-zero source point, placed back onto the NxN frequency grid.

2. **`read_mask_oas(filepath, sys)`** — reads OASIS via `gdstk`. Uses `matplotlib.path.Path.contains_points` to rasterize each polygon onto pixel centres `((col+0.5)*dx, (row+0.5)*dx)` — supports arbitrary polygon shapes including rotated layouts. Bounding box is used only to restrict which pixels to test, not to fill.

3. **`make_grids` / `make_pupil`** — identical to `1_basic`. Pupil includes a fixed Zernike coma term `(n=3, m=1): 0.02 waves`.

4. **`abbe_simulation` / `hopkins_simulation`** — identical to `1_basic`. `--n_svd` controls Hopkins SVD truncation (default 30).

5. **Plotting** — `plot_results` (6-panel summary) and `plot_linecuts` (horizontal linecut at a configurable y position), saved as `<PREFIX>_aerial_image.png` and `<PREFIX>_linecut.png`. `plot_linecuts` accepts `y_nm` (default 0 = centre row); the nearest grid row is selected and the actual y is shown in the plot title.

### `generate_layout.py` — layout generator

Creates an OASIS file with 16 L-shaped polygons (arm width=32 nm, arm length=96 nm) at random orientations on a 4×4 grid (pitch=220 nm, origin offset=150 nm). The 150 nm offset ensures rotated L-shapes (max radius ~68 nm) stay fully within the 1024×1024 nm simulation field.

Key functions: `l_shape_vertices(w, h)` returns the 6 canonical vertices centred at origin; `rotate(pts, angle_deg)` applies a 2D rotation matrix; `translate(pts, dx, dy)` shifts to the grid position.

### Unit convention

All lengths in **nm**, all spatial frequencies in **1/nm** — consistent with `1_basic`.

### Relationship to 1_basic

`1_basic/aerial_image_sim.py` generates masks and source in-code then exports them. This folder consumes those exports. `source_distribution.csv` is copied from `1_basic` and bundled here for convenience.

## Example output: multi_L-shaped

16 L-shapes (arm width=32 nm, arm length=96 nm) at random orientations, placed on a 4×4 grid (pitch=220 nm, offset=150 nm), filling the 1024×1024 nm field uniformly. Simulated with λ=13.5 nm EUV, NA=0.33, annular source (σ=0.55–0.8), 256×256 grid, dx=4 nm.

**`multi_L-shaped_aerial_image.png`** — 6-panel summary showing the rasterized mask with 16 distinct L orientations and the corresponding Abbe and Hopkins aerial images. Each L appears as a diffraction-blurred spot with an asymmetric L profile. Mask transmission is ~7.8%. Abbe–Hopkins normalised RMS = 0.00097.

![](multi_L-shaped_aerial_image.png)

**`multi_L-shaped_linecut.png`** — horizontal cut at y=100 nm, which passes through row 3 of L-shapes (centred at y≈78 nm OAS, just 22 nm below the cut). Four well-separated intensity peaks appear at x≈−362, −142, +78, +298 nm, corresponding to the 4 columns. Peak heights vary (0.8–1.0 normalised) because each L has a different orientation — orientations that expose more chrome along this row produce taller peaks. The rightmost peak shows a shoulder/double structure from an L whose both arms contribute at this y position. Abbe and Hopkins are nearly indistinguishable throughout.

![](multi_L-shaped_linecut.png)
