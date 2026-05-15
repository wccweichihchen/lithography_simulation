# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the simulation

```bash
python aerial_image_sim.py
```

This runs all four mask patterns sequentially and writes outputs to the current directory. Requires `numpy`, `matplotlib`, and `gdstk` (or `gdspy` as fallback for OASIS export).

## Architecture

The entire simulation lives in `aerial_image_sim.py`, structured as a pipeline of independent modules:

1. **`OpticalSystem` dataclass** — single source of truth for all optical parameters (wavelength, NA, sigma, grid size). Passed by reference through every stage. All lengths are in **nanometres**; all spatial frequencies in **1/nm**.

2. **Coordinate grids** (`make_grids`) — produces spatial `(X, Y)` and frequency `(FX, FY)` meshgrids, both centered (fftshift convention). Used by source, pupil, mask, and both simulation methods.

3. **Source** (`make_source`) — discrete illumination pupil sampled on the frequency grid. Supports `conventional`, `annular`, `dipole`, and `quasar` shapes. Normalised so `sum == 1`.

4. **Pupil** (`make_pupil`) — coherent transfer function (complex-valued). Supports Zernike aberrations via `zernike_coeffs={(n, m): coeff_waves}`.

5. **Mask** (`make_mask`) — binary amplitude mask, four patterns: `lines`, `contact`, `L-shaped`, `checkerboard`. CD and pitch in nm.

6. **Simulation methods** — two independent implementations of partially coherent aerial image:
   - **Abbe** (`abbe_simulation`): direct incoherent sum over source points. Exact but O(N_source) FFTs.
   - **Hopkins** (`hopkins_simulation`): SVD decomposition of the Transmission Cross Coefficient (TCC). `n_svd` controls truncation; default 30 terms. More efficient for many source points.

7. **Outputs** — per-pattern PNG aerial images and linecut plots, `.oas` mask files (OASIS via gdstk), and `source_distribution.csv` (sparse source point export).

## Key conventions

- The fftshift/ifftshift pair wraps every FFT/IFFT call; the pupil and source are stored in shifted (DC-centered) form.
- Abbe shifts the mask spectrum by rolling the array; Hopkins shifts the pupil instead (sign is negated).
- `write_oasis` uses run-length + row-merge compression to minimise rectangle count before writing.
- If `gdstk` is unavailable, `write_oasis` falls back to `gdspy` and writes GDSII (`.gds`) instead of OASIS (`.oas`).
- If `/mnt/user-data/outputs` exists (container environment), all output files are copied there at the end.
