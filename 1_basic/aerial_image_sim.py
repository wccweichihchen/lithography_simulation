"""
Lithography Aerial Image Simulation
=====================================
Implements both Abbe and Hopkins methods for aerial image computation.

Physical model:
  - Partially coherent illumination (scalar approximation)
  - Thin mask approximation
  - Circular pupil with optional Zernike aberrations
  - Abbe method  : incoherent sum over source points
  - Hopkins method: Transmission Cross Coefficient (TCC) via SVD

Unit convention (throughout the entire code):
  - All lengths            : nanometres [nm]
  - All spatial frequencies: [1/nm]

Usage:
  python aerial_image_sim.py
"""

import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from dataclasses import dataclass
from typing import Literal

# ---------------------------------------------------------------------------
# 1.  Optical system parameters
# ---------------------------------------------------------------------------

@dataclass
class OpticalSystem:
    """Scalar partially-coherent optical system.  All lengths in nm."""
    wavelength  : float = 13.5    # [nm]  EUV
    NA          : float = 0.33    # numerical aperture
    sigma       : float = 0.8     # partial coherence (outer)
    sigma_inner : float = 0.55    # inner sigma (0 = conventional)
    source_type : Literal["conventional", "annular", "dipole", "quasar"] = "annular"
    N           : int   = 256     # grid points per side
    dx          : float = 4.0     # [nm] spatial sampling

    @property
    def f_max(self): return self.NA / self.wavelength   # [1/nm]
    @property
    def df(self):    return 1.0 / (self.N * self.dx)   # [1/nm]


# ---------------------------------------------------------------------------
# 2.  Coordinate grids
# ---------------------------------------------------------------------------

def make_grids(sys: OpticalSystem):
    """Return spatial [nm] and frequency [1/nm] coordinate grids."""
    N, dx = sys.N, sys.dx
    x  = (np.arange(N) - N // 2) * dx               # [nm]
    fx = np.fft.fftshift(np.fft.fftfreq(N, d=dx))   # [1/nm]
    X,  Y  = np.meshgrid(x,  x,  indexing='xy')
    FX, FY = np.meshgrid(fx, fx, indexing='xy')
    return X, Y, FX, FY


# ---------------------------------------------------------------------------
# 3.  Light source
# ---------------------------------------------------------------------------

def make_source(sys: OpticalSystem, FX, FY):
    """
    Discrete source distribution J(fx, fy), normalised so sum == 1.
    rho = sqrt(fx^2 + fy^2) / (NA/lambda)  ->  rho=1 at outer sigma edge.
    """
    norm = sys.NA / sys.wavelength   # [1/nm]
    rho  = np.sqrt(FX**2 + FY**2) / norm
    phi  = np.arctan2(FY, FX)
    src  = np.zeros_like(rho)

    if sys.source_type == "conventional":
        src[rho <= sys.sigma] = 1.0

    elif sys.source_type == "annular":
        src[(rho >= sys.sigma_inner) & (rho <= sys.sigma)] = 1.0

    elif sys.source_type == "dipole":
        spot_half   = 0.15
        pole_offset = (sys.sigma + max(sys.sigma_inner, 0.1)) / 2
        for angle in [0, np.pi]:
            cx = pole_offset * np.cos(angle) / sys.sigma
            cy = pole_offset * np.sin(angle) / sys.sigma
            d  = np.sqrt((rho * np.cos(phi) - cx)**2 +
                         (rho * np.sin(phi) - cy)**2)
            src[d <= spot_half] = 1.0

    elif sys.source_type == "quasar":
        spot_half   = 0.15
        pole_offset = (sys.sigma + max(sys.sigma_inner, 0.1)) / 2
        for angle in [np.pi/4, 3*np.pi/4, 5*np.pi/4, 7*np.pi/4]:
            cx = pole_offset * np.cos(angle) / sys.sigma
            cy = pole_offset * np.sin(angle) / sys.sigma
            d  = np.sqrt((rho * np.cos(phi) - cx)**2 +
                         (rho * np.sin(phi) - cy)**2)
            src[d <= spot_half] = 1.0

    src[rho > sys.sigma] = 0.0
    total = src.sum()
    if total > 0:
        src /= total
    return src


# ---------------------------------------------------------------------------
# 4.  Pupil function
# ---------------------------------------------------------------------------

def make_pupil(sys: OpticalSystem, FX, FY, zernike_coeffs=None):
    """
    Coherent transfer function (pupil).
    zernike_coeffs: dict {(n, m): coefficient [waves]}
    """
    norm     = sys.NA / sys.wavelength
    rho_lens = np.sqrt(FX**2 + FY**2) / norm
    phi_lens = np.arctan2(FY, FX)
    in_pupil = rho_lens <= 1.0
    phase    = np.zeros_like(rho_lens)

    if zernike_coeffs:
        for (n, m), coeff in zernike_coeffs.items():
            phase += coeff * _zernike(n, m, rho_lens, phi_lens)

    H = np.zeros_like(rho_lens, dtype=complex)
    H[in_pupil] = np.exp(1j * 2 * np.pi * phase[in_pupil])
    return H


def _zernike(n, m, rho, phi):
    R = _radial_zernike(n, abs(m), rho)
    if   m == 0: return R
    elif m  > 0: return R * np.cos(m * phi)
    else:        return R * np.sin(abs(m) * phi)


def _radial_zernike(n, m, rho):
    R = np.zeros_like(rho)
    for s in range((n - m) // 2 + 1):
        c = ((-1)**s * math.factorial(n - s) /
             (math.factorial(s) *
              math.factorial((n + m) // 2 - s) *
              math.factorial((n - m) // 2 - s)))
        R += c * rho**(n - 2 * s)
    return R


# ---------------------------------------------------------------------------
# 5.  Mask patterns
# ---------------------------------------------------------------------------

def make_mask(sys: OpticalSystem, X, Y,
              pattern: Literal["lines", "contact", "L-shaped", "checkerboard"] = "lines",
              cd: float = 32.0, pitch: float = 64.0):
    """
    Binary chrome-on-glass mask (amplitude: 0 = clear, 1 = chrome).

    Parameters
    ----------
    cd    : critical dimension [nm]  - line / hole width
    pitch : period [nm]
    X, Y  : coordinate arrays [nm]
    """
    mask = np.zeros((sys.N, sys.N), dtype=float)

    if pattern == "lines":
        # Vertical lines: chrome width = cd, period = pitch
        mask[np.mod(X, pitch) < cd] = 1.0

    elif pattern == "contact":
        # 2-D periodic square contact holes: cd x cd on pitch x pitch grid
        mask[(np.mod(X, pitch) < cd) & (np.mod(Y, pitch) < cd)] = 1.0

    elif pattern == "L-shaped":
        # Single L-shape centred in field.
        # Arm width = cd, arm length = 3*cd. L opens toward top-right.
        cx  = sys.N // 2
        cy  = sys.N // 2
        w   = int(cd / sys.dx)      # arm width in pixels
        arm = 3                     # arm length = arm * w pixels
        x0  = cx - arm * w // 2    # bounding-box origin
        y0  = cy - arm * w // 2
        x1  = x0 + arm * w
        y1  = y0 + arm * w
        mask[y0:y1,   x0:x0+w] = 1.0   # vertical arm (left side)
        mask[y0:y0+w, x0:x1  ] = 1.0   # horizontal arm (bottom)

    elif pattern == "checkerboard":
        # Alternating squares; each square is pitch x pitch
        x_q = (np.mod(X, 2 * pitch) < pitch)
        y_q = (np.mod(Y, 2 * pitch) < pitch)
        mask[x_q == y_q] = 1.0

    return mask


# ---------------------------------------------------------------------------
# 6.  Abbe method
# ---------------------------------------------------------------------------

def abbe_simulation(sys: OpticalSystem, mask, source, pupil, FX, FY):
    """
    Abbe aerial image: incoherent sum over source points.

      I(x) = sum_s  J(s) |IFFT[ M_tilde(f - s) * H(f) ]|^2
    """
    N    = sys.N
    nz   = np.nonzero(source.ravel())[0]
    print(f"    Abbe: {len(nz)} source points")

    M_fft = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(mask)))
    AI    = np.zeros((N, N), dtype=float)

    for idx in nz:
        j_val = source.ravel()[idx]
        dix   = int(round(FX.ravel()[idx] / sys.df))
        diy   = int(round(FY.ravel()[idx] / sys.df))
        M_sh  = np.roll(np.roll(M_fft, dix, axis=1), diy, axis=0)
        field = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(M_sh * pupil)))
        AI   += j_val * np.abs(field)**2

    return AI


# ---------------------------------------------------------------------------
# 7.  Hopkins method (TCC via SVD)
# ---------------------------------------------------------------------------

def hopkins_simulation(sys: OpticalSystem, mask, source, pupil, FX, FY,
                       n_svd: int = 30):
    """
    Hopkins aerial image via SVD decomposition of the TCC.

    TCC(f1, f2) = sum_s J(s) H(s+f1) H*(s+f2)
    SVD: TCC ~= sum_k lambda_k phi_k(f1) phi_k*(f2)
    AI(x)      = sum_k lambda_k |IFFT[M_tilde(f) * phi_k(f)]|^2
    """
    N    = sys.N
    nz   = np.nonzero(source.ravel())[0]
    Ns   = len(nz)
    print(f"    Hopkins: {Ns} source pts, {n_svd} SVD terms")

    # Build A using ALL source points so no energy is lost.
    # A[:,s] = sqrt(J(s)) * vec(pupil shifted by source point s)
    # TCC = A A†  ->  SVD of A gives eigenvectors of TCC.
    # Truncation to n_svd terms happens AFTER the SVD, not before.
    A = np.zeros((N * N, Ns), dtype=complex)

    for col, idx in enumerate(nz):
        j_val = source.ravel()[idx]
        dix   = int(round(FX.ravel()[idx] / sys.df))
        diy   = int(round(FY.ravel()[idx] / sys.df))
        H_sh  = np.roll(np.roll(pupil, -dix, axis=1), -diy, axis=0)
        A[:, col] = np.sqrt(j_val) * H_sh.ravel()

    U, s_vals, _ = np.linalg.svd(A, full_matrices=False)
    eigenvalues  = s_vals**2
    eigenvectors = U

    # Truncate to n_svd dominant terms AFTER full SVD
    r = min(n_svd, len(eigenvalues))
    eigenvalues  = eigenvalues[:r]
    eigenvectors = eigenvectors[:, :r]

    M_fft = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(mask)))
    AI    = np.zeros((N, N), dtype=float)

    for k in range(r):
        phi_k  = eigenvectors[:, k].reshape(N, N)
        kernel = M_fft * phi_k
        img_k  = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(kernel)))
        AI    += eigenvalues[k] * np.abs(img_k)**2

    return AI


# ---------------------------------------------------------------------------
# 8.  Plotting
# ---------------------------------------------------------------------------

def plot_results(sys: OpticalSystem, mask, source, pupil,
                 AI_abbe, AI_hop, X, Y, FX, FY, pattern_name: str = ""):
    """Full summary figure: mask, source, pupil, both aerial images."""
    fig = plt.figure(figsize=(18, 11))
    fig.patch.set_facecolor('#0f0f1a')
    gs  = gridspec.GridSpec(2, 4, figure=fig, hspace=0.38, wspace=0.35)

    norm_f = sys.NA / sys.wavelength   # [1/nm]
    ext_x  = [X[0, 0],  X[0, -1],  Y[0, 0],  Y[-1, 0]]   # [nm]
    ext_f  = [FX[0, 0] / norm_f, FX[0, -1] / norm_f,
              FY[0, 0] / norm_f, FY[-1, 0] / norm_f]

    def styled(ax, title):
        ax.set_facecolor('#0f0f1a')
        ax.set_title(title, color='white', fontsize=10, pad=6)
        for sp in ax.spines.values(): sp.set_color('#444')
        ax.tick_params(colors='#aaa', labelsize=7)

    def cbar(im, ax):
        cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.ax.yaxis.set_tick_params(color='#aaa', labelsize=7)
        plt.setp(cb.ax.yaxis.get_ticklabels(), color='#aaa')

    # Mask
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(mask, extent=ext_x, cmap='inferno', origin='lower')
    cbar(im, ax)
    ax.set_xlabel('x [nm]', color='#aaa', fontsize=8)
    ax.set_ylabel('y [nm]', color='#aaa', fontsize=8)
    styled(ax, f'Mask - {pattern_name}')

    # Source
    ax = fig.add_subplot(gs[0, 1])
    vmin = max(source.max() * 1e-4, 1e-20)
    im = ax.imshow(source, extent=ext_f, cmap='hot', origin='lower',
                   norm=LogNorm(vmin=vmin, vmax=source.max()))
    cbar(im, ax)
    ax.set_xlabel('sigma_x', color='#aaa', fontsize=8)
    ax.set_ylabel('sigma_y', color='#aaa', fontsize=8)
    styled(ax, f'Source [{sys.source_type}]  sigma={sys.sigma}')

    # Pupil magnitude
    ax = fig.add_subplot(gs[0, 2])
    im = ax.imshow(np.abs(pupil), extent=ext_f, cmap='viridis', origin='lower')
    cbar(im, ax)
    ax.set_xlabel('fx / (NA/lam)', color='#aaa', fontsize=8)
    styled(ax, f'|Pupil|  NA={sys.NA}')

    # Pupil phase
    ax = fig.add_subplot(gs[0, 3])
    ph = np.angle(pupil) * (np.abs(pupil) > 0.01)
    im = ax.imshow(ph, extent=ext_f, cmap='hsv', origin='lower',
                   vmin=-np.pi, vmax=np.pi)
    cbar(im, ax)
    ax.set_xlabel('fx / (NA/lam)', color='#aaa', fontsize=8)
    styled(ax, 'Pupil phase [rad]')

    # Abbe aerial image
    ax = fig.add_subplot(gs[1, 0:2])
    im = ax.imshow(AI_abbe, extent=ext_x, cmap='inferno', origin='lower')
    cbar(im, ax)
    ax.set_xlabel('x [nm]', color='#aaa', fontsize=8)
    ax.set_ylabel('y [nm]', color='#aaa', fontsize=8)
    styled(ax, 'Aerial Image - Abbe')

    # Hopkins aerial image
    ax = fig.add_subplot(gs[1, 2:4])
    im = ax.imshow(AI_hop, extent=ext_x, cmap='inferno', origin='lower')
    cbar(im, ax)
    ax.set_xlabel('x [nm]', color='#aaa', fontsize=8)
    styled(ax, 'Aerial Image - Hopkins (SVD TCC)')

    fig.suptitle(
        f'Aerial Image Simulation - {pattern_name}\n'
        f'lam={sys.wavelength:.1f} nm  NA={sys.NA}  sigma={sys.sigma}'
        f'  {sys.N}x{sys.N}  dx={sys.dx:.1f} nm',
        color='white', fontsize=12, y=1.01)

    fname = f'aerial_image_{pattern_name}.png'
    plt.savefig(fname, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    print(f"    Saved: {fname}")
    plt.close()


def plot_linecuts(AI_abbe, AI_hop, sys: OpticalSystem, pattern_name: str = ""):
    """Horizontal line-cut through image centre."""
    N  = sys.N
    cx = N // 2
    x  = (np.arange(N) - cx) * sys.dx   # [nm]

    a_cut = AI_abbe[cx, :] / AI_abbe.max()
    h_cut = AI_hop [cx, :] / AI_hop.max()

    fig, ax = plt.subplots(figsize=(9, 4))
    fig.patch.set_facecolor('#0f0f1a')
    ax.set_facecolor('#0f0f1a')
    ax.plot(x, a_cut, color='#ff6b35', lw=2,        label='Abbe')
    ax.plot(x, h_cut, color='#4ecdc4', lw=2, ls='--', label='Hopkins')
    ax.set_xlabel('x [nm]', color='#aaa')
    ax.set_ylabel('Normalised intensity', color='#aaa')
    ax.set_title(f'Line-cut (centre row) - {pattern_name}',
                 color='white', fontsize=11)
    ax.legend(facecolor='#1a1a2e', labelcolor='white')
    ax.tick_params(colors='#aaa')
    for sp in ax.spines.values(): sp.set_color('#444')
    ax.grid(color='#333', lw=0.5)

    plt.tight_layout()
    fname = f'aerial_linecut_{pattern_name}.png'
    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='#0f0f1a')
    print(f"    Saved: {fname}")
    plt.close()


# ---------------------------------------------------------------------------
# 9.  OASIS mask export  (uses gdstk)
# ---------------------------------------------------------------------------

def write_oasis(mask: np.ndarray, sys: OpticalSystem, filepath: str,
                layer: int = 1, datatype: int = 0):
    """
    Write binary mask to OASIS using gdstk (pip install gdstk).
    Falls back to gdspy -> GDSII if gdstk is unavailable.

    Coordinate convention
    ---------------------
    Units  : nm  (1 gdstk user unit = 1 nm)
    Origin : (0, 0) at bottom-left corner of the mask array
    row 0  -> y = 0          (bottom of layout)
    row N-1-> y = (N-1)*dx   (top of layout)
    """
    try:
        import gdstk
        _backend = 'gdstk'
    except ImportError:
        try:
            import gdspy as gdstk
            _backend = 'gdspy'
        except ImportError:
            raise ImportError(
                "Neither gdstk nor gdspy is installed.\n"
                "Install with:  pip install gdstk")

    dx  = sys.dx   # [nm]
    bm  = (mask >= 0.5)
    N   = bm.shape[0]

    # ---- Run-length encode each row ----------------------------------------
    def row_runs(row):
        runs = []; in_run = False
        for x, v in enumerate(row):
            if v and not in_run:
                xs = x; in_run = True
            elif not v and in_run:
                runs.append((xs, x)); in_run = False
        if in_run:
            runs.append((xs, len(row)))
        return frozenset(runs)

    all_runs = [row_runs(bm[row]) for row in range(N)]

    # ---- Merge consecutive identical rows -> rectangles --------------------
    active = {}   # seg -> y_bot when this run opened
    rects  = []   # (x0, y_bot, x1, y_top) in nm

    for row_idx in range(N):
        y_bot = row_idx       * dx
        y_top = (row_idx + 1) * dx
        cur   = all_runs[row_idx]
        prev  = set(active.keys())

        for seg in prev - cur:                        # runs that just closed
            rects.append((seg[0] * dx, active.pop(seg),
                          seg[1] * dx, y_bot))

        for seg in cur - prev:                        # runs that just opened
            active[seg] = y_bot

    for seg, y_open in active.items():                # flush still-open runs
        rects.append((seg[0] * dx, y_open, seg[1] * dx, N * dx))

    print(f"    OASIS: {len(rects)} rectangles  [{_backend}]")

    # ---- Write via gdstk ---------------------------------------------------
    lib  = gdstk.Library() if _backend == 'gdstk' else gdstk.GdsLibrary()
    cell = lib.new_cell('TOP')

    for (x0, y0, x1, y1) in rects:
        rect = gdstk.rectangle((x0, y0), (x1, y1),
                               layer=layer, datatype=datatype)
        cell.add(rect)

    if _backend == 'gdstk':
        lib.write_oas(filepath)
    else:
        gds_path = filepath.replace('.oas', '.gds')
        lib.write_gds(gds_path)
        filepath = gds_path
        print("    Note: gdspy has no OASIS writer; wrote GDSII instead.")

    import os
    print(f"    Saved: {filepath}  ({os.path.getsize(filepath)/1024:.1f} kB)")


# ---------------------------------------------------------------------------
# 10. Source CSV export
# ---------------------------------------------------------------------------

def write_source_csv(source: np.ndarray, sys: OpticalSystem, filepath: str):
    """
    Export non-zero source points to CSV (sparse format).

    Columns
    -------
    sigma_x, sigma_y : normalised pupil coords  (1 = outer sigma edge)
    rho, phi_deg     : polar form
    fx_per_nm        : spatial frequency x  [1/nm]
    fy_per_nm        : spatial frequency y  [1/nm]
    intensity        : weight (sum over all points = 1)
    """
    import csv

    norm  = sys.NA / sys.wavelength
    N     = sys.N
    fx_ax = (np.arange(N) - N // 2) * sys.df   # [1/nm]

    rows = []
    for iy in range(N):
        for ix in range(N):
            val = source[iy, ix]
            if val == 0.0:
                continue
            fx  = fx_ax[ix];  fy = fx_ax[iy]
            sx  = fx / norm;  sy = fy / norm
            rho = np.sqrt(sx**2 + sy**2)
            phi = np.degrees(np.arctan2(sy, sx))
            rows.append({
                'sigma_x'   : f'{sx:.6f}',
                'sigma_y'   : f'{sy:.6f}',
                'rho'       : f'{rho:.6f}',
                'phi_deg'   : f'{phi:.4f}',
                'fx_per_nm' : f'{fx:.6e}',
                'fy_per_nm' : f'{fy:.6e}',
                'intensity' : f'{val:.8e}',
            })

    header = [
        '# Lithography Source Distribution',
        f'# wavelength_nm = {sys.wavelength:.3f}',
        f'# NA            = {sys.NA}',
        f'# sigma_outer   = {sys.sigma}',
        f'# sigma_inner   = {sys.sigma_inner}',
        f'# source_type   = {sys.source_type}',
        f'# grid_N        = {N}',
        f'# dx_nm         = {sys.dx:.3f}',
        f'# df_per_nm     = {sys.df:.6e}',
        f'# total_points  = {len(rows)}',
        f'# norm_conv     : sigma_x = fx_per_nm / (NA/wavelength)',
        f'# intensity_sum = {source.sum():.6f}  (should be 1.0)',
    ]

    with open(filepath, 'w', newline='\n') as f:
        for line in header:
            f.write(line + '\n')
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

    print(f"    Saved: {filepath}  ({len(rows)} source points)")


# ---------------------------------------------------------------------------
# 11.  Main - run all four mask types
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import os, shutil

    print("=" * 60)
    print("  Lithography Aerial Image Simulation  (units: nm)")
    print("=" * 60)

    # -- Optical system (all lengths in nm) ----------------------------------
    sys = OpticalSystem(
        wavelength  = 13.5,     # [nm]  EUV
        NA          = 0.33,
        sigma       = 0.8,
        sigma_inner = 0.55,
        source_type = "annular",
        N           = 256,
        dx          = 4.0,      # [nm]
    )

    print(f"\n  Optical params:")
    print(f"    wavelength = {sys.wavelength:.1f} nm")
    print(f"    NA         = {sys.NA}")
    print(f"    sigma      = {sys.sigma}  (inner = {sys.sigma_inner})")
    print(f"    source     = {sys.source_type}")
    print(f"    grid       = {sys.N}x{sys.N},  dx = {sys.dx:.1f} nm")
    print(f"    FoV        = {sys.N * sys.dx:.0f} nm")

    # -- Shared grids, source, pupil -----------------------------------------
    X, Y, FX, FY = make_grids(sys)

    print("\n  Building source & pupil...")
    source = make_source(sys, FX, FY)
    print(f"    Non-zero source points: {np.count_nonzero(source)}")
    pupil  = make_pupil(sys, FX, FY, zernike_coeffs={(3, 1): 0.02})

    print("\n  Exporting source to CSV...")
    write_source_csv(source, sys, 'source_distribution.csv')

    # -- Mask configurations -------------------------------------------------
    mask_configs = [
        dict(pattern="lines",        cd=32.0, pitch=64.0),
        dict(pattern="contact",      cd=32.0, pitch=64.0),
        dict(pattern="checkerboard", cd=32.0, pitch=64.0),
        dict(pattern="L-shaped",     cd=32.0, pitch=64.0),
    ]

    output_files = ['source_distribution.csv']

    # -- Loop over all mask types --------------------------------------------
    for cfg in mask_configs:
        pattern = cfg["pattern"]
        print(f"\n{'─'*60}")
        print(f"  Pattern : {pattern}  "
              f"cd={cfg['cd']:.0f} nm  pitch={cfg['pitch']:.0f} nm")
        print(f"{'─'*60}")

        mask = make_mask(sys, X, Y, **cfg)
        print(f"  Mask transmission : {mask.mean():.3f}")

        print("  Running Abbe simulation...")
        AI_abbe = abbe_simulation(sys, mask, source, pupil, FX, FY)
        print(f"    Intensity range : [{AI_abbe.min():.4f}, {AI_abbe.max():.4f}]")

        print("  Running Hopkins simulation...")
        AI_hop = hopkins_simulation(sys, mask, source, pupil, FX, FY, n_svd=30)
        print(f"    Intensity range : [{AI_hop.min():.4f}, {AI_hop.max():.4f}]")

        a_n = AI_abbe / AI_abbe.max()
        h_n = AI_hop  / AI_hop.max()
        print(f"  Norm. RMS(Abbe-Hopkins) = {np.sqrt(np.mean((a_n - h_n)**2)):.5f}")

        print("  Plotting...")
        plot_results(sys, mask, source, pupil, AI_abbe, AI_hop,
                     X, Y, FX, FY, pattern_name=pattern)
        plot_linecuts(AI_abbe, AI_hop, sys, pattern_name=pattern)

        print("  Exporting OASIS...")
        oas_file = f'mask_{pattern}.oas'
        try:
            write_oasis(mask, sys, filepath=oas_file, layer=1, datatype=0)
            output_files.append(oas_file)
        except ImportError as e:
            print(f"    Skipped (install gdstk): {e}")

        output_files += [
            f'aerial_image_{pattern}.png',
            f'aerial_linecut_{pattern}.png',
        ]

    # -- Copy all outputs to /mnt/user-data/outputs (container env) ----------
    out_dir = '/mnt/user-data/outputs'
    if os.path.isdir(out_dir):
        for fname in output_files:
            if os.path.exists(fname):
                shutil.copy(fname, os.path.join(out_dir, fname))
                print(f"  Copied -> {out_dir}/{fname}")

    print(f"\n{'='*60}")
    print("  All patterns complete.")
    print(f"{'='*60}")
