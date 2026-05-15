"""
Lithography Aerial Image Simulation — from source CSV + mask OASIS
==================================================================
Reads a pre-computed source distribution (CSV) and a binary mask (OASIS),
then computes aerial images using both Abbe and Hopkins methods.

Usage:
    python sim.py --source source.csv --mask mask.oas [--n_svd 30] [--output PREFIX]

All lengths: nm.  All spatial frequencies: 1/nm.
"""

import argparse
import csv
import io
import math
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Optical system
# ---------------------------------------------------------------------------

@dataclass
class OpticalSystem:
    wavelength  : float = 13.5
    NA          : float = 0.33
    sigma       : float = 0.8
    sigma_inner : float = 0.55
    source_type : str   = "annular"
    N           : int   = 256
    dx          : float = 4.0

    @property
    def f_max(self): return self.NA / self.wavelength
    @property
    def df(self):    return 1.0 / (self.N * self.dx)


# ---------------------------------------------------------------------------
# Coordinate grids
# ---------------------------------------------------------------------------

def make_grids(sys):
    N, dx = sys.N, sys.dx
    x  = (np.arange(N) - N // 2) * dx
    fx = np.fft.fftshift(np.fft.fftfreq(N, d=dx))
    X,  Y  = np.meshgrid(x,  x,  indexing='xy')
    FX, FY = np.meshgrid(fx, fx, indexing='xy')
    return X, Y, FX, FY


# ---------------------------------------------------------------------------
# Pupil
# ---------------------------------------------------------------------------

def make_pupil(sys, FX, FY, zernike_coeffs=None):
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
# Source reader
# ---------------------------------------------------------------------------

def read_source_csv(filepath):
    """
    Read source CSV written by write_source_csv() in 1_basic.
    Returns (source_array [N x N], OpticalSystem).

    The CSV has comment lines starting with '#' carrying optical metadata,
    followed by a header row and data rows with columns:
        sigma_x, sigma_y, rho, phi_deg, fx_per_nm, fy_per_nm, intensity
    """
    meta = {}
    data_lines = []

    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                parts = line.lstrip('#').strip().split('=')
                if len(parts) == 2:
                    meta[parts[0].strip()] = parts[1].strip()
            else:
                data_lines.append(line)

    reader = csv.DictReader(io.StringIO(''.join(data_lines)))
    rows = list(reader)

    wavelength  = float(meta.get('wavelength_nm', 13.5))
    NA          = float(meta.get('NA', 0.33))
    sigma       = float(meta.get('sigma_outer', 0.8))
    sigma_inner = float(meta.get('sigma_inner', 0.55))
    source_type = meta.get('source_type', 'annular')
    N           = int(meta.get('grid_N', 256))
    dx          = float(meta.get('dx_nm', 4.0))

    sys = OpticalSystem(
        wavelength=wavelength, NA=NA,
        sigma=sigma, sigma_inner=sigma_inner,
        source_type=source_type, N=N, dx=dx,
    )

    df     = 1.0 / (N * dx)
    source = np.zeros((N, N), dtype=float)

    for row in rows:
        fx  = float(row['fx_per_nm'])
        fy  = float(row['fy_per_nm'])
        val = float(row['intensity'])
        ix  = int(round(fx / df)) + N // 2
        iy  = int(round(fy / df)) + N // 2
        if 0 <= ix < N and 0 <= iy < N:
            source[iy, ix] = val

    return source, sys


# ---------------------------------------------------------------------------
# Mask reader
# ---------------------------------------------------------------------------

def read_mask_oas(filepath, sys):
    """
    Read OASIS mask written by write_oasis() in 1_basic.
    Returns binary mask array of shape (N, N).

    Coordinate convention (matches write_oasis):
        OAS x in [0, N*dx], OAS y in [0, N*dx]
        col = round(x / dx), row = round(y / dx)
    """
    try:
        import gdstk
    except ImportError:
        raise ImportError("gdstk is required to read OASIS files: pip install gdstk")

    lib  = gdstk.read_oas(filepath)
    tops = lib.top_level()
    if not tops:
        raise ValueError(f"No top-level cell in {filepath}")
    cell = tops[0]

    N  = sys.N
    dx = sys.dx
    mask = np.zeros((N, N), dtype=float)

    for poly in cell.polygons:
        pts = poly.points
        x0  = pts[:, 0].min()
        x1  = pts[:, 0].max()
        y0  = pts[:, 1].min()
        y1  = pts[:, 1].max()
        c0  = max(0, int(round(x0 / dx)))
        c1  = min(N, int(round(x1 / dx)))
        r0  = max(0, int(round(y0 / dx)))
        r1  = min(N, int(round(y1 / dx)))
        mask[r0:r1, c0:c1] = 1.0

    return mask


# ---------------------------------------------------------------------------
# Abbe simulation
# ---------------------------------------------------------------------------

def abbe_simulation(sys, mask, source, pupil, FX, FY):
    N    = sys.N
    nz   = np.nonzero(source.ravel())[0]
    print(f"  Abbe: {len(nz)} source points")

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
# Hopkins simulation (TCC via SVD)
# ---------------------------------------------------------------------------

def hopkins_simulation(sys, mask, source, pupil, FX, FY, n_svd=30):
    N    = sys.N
    nz   = np.nonzero(source.ravel())[0]
    Ns   = len(nz)
    print(f"  Hopkins: {Ns} source pts, {n_svd} SVD terms")

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

    r = min(n_svd, len(eigenvalues))
    eigenvalues  = eigenvalues[:r]
    eigenvectors = eigenvectors[:, :r]

    M_fft = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(mask)))
    AI    = np.zeros((N, N), dtype=float)

    for k in range(r):
        phi_k  = eigenvectors[:, k].reshape(N, N)
        img_k  = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(M_fft * phi_k)))
        AI    += eigenvalues[k] * np.abs(img_k)**2

    return AI


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_results(sys, mask, source, pupil, AI_abbe, AI_hop,
                 X, Y, FX, FY, output_prefix):
    fig = plt.figure(figsize=(18, 11))
    fig.patch.set_facecolor('#0f0f1a')
    gs  = gridspec.GridSpec(2, 4, figure=fig, hspace=0.38, wspace=0.35)

    norm_f = sys.NA / sys.wavelength
    ext_x  = [X[0, 0], X[0, -1], Y[0, 0], Y[-1, 0]]
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

    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(mask, extent=ext_x, cmap='inferno', origin='lower')
    cbar(im, ax)
    ax.set_xlabel('x [nm]', color='#aaa', fontsize=8)
    ax.set_ylabel('y [nm]', color='#aaa', fontsize=8)
    styled(ax, 'Mask')

    ax = fig.add_subplot(gs[0, 1])
    vmin = max(source.max() * 1e-4, 1e-20)
    im = ax.imshow(source, extent=ext_f, cmap='hot', origin='lower',
                   norm=LogNorm(vmin=vmin, vmax=source.max()))
    cbar(im, ax)
    ax.set_xlabel('sigma_x', color='#aaa', fontsize=8)
    ax.set_ylabel('sigma_y', color='#aaa', fontsize=8)
    styled(ax, f'Source [{sys.source_type}]  sigma={sys.sigma}')

    ax = fig.add_subplot(gs[0, 2])
    im = ax.imshow(np.abs(pupil), extent=ext_f, cmap='viridis', origin='lower')
    cbar(im, ax)
    ax.set_xlabel('fx / (NA/lam)', color='#aaa', fontsize=8)
    styled(ax, f'|Pupil|  NA={sys.NA}')

    ax = fig.add_subplot(gs[0, 3])
    ph = np.angle(pupil) * (np.abs(pupil) > 0.01)
    im = ax.imshow(ph, extent=ext_f, cmap='hsv', origin='lower',
                   vmin=-np.pi, vmax=np.pi)
    cbar(im, ax)
    ax.set_xlabel('fx / (NA/lam)', color='#aaa', fontsize=8)
    styled(ax, 'Pupil phase [rad]')

    ax = fig.add_subplot(gs[1, 0:2])
    im = ax.imshow(AI_abbe, extent=ext_x, cmap='inferno', origin='lower')
    cbar(im, ax)
    ax.set_xlabel('x [nm]', color='#aaa', fontsize=8)
    ax.set_ylabel('y [nm]', color='#aaa', fontsize=8)
    styled(ax, 'Aerial Image - Abbe')

    ax = fig.add_subplot(gs[1, 2:4])
    im = ax.imshow(AI_hop, extent=ext_x, cmap='inferno', origin='lower')
    cbar(im, ax)
    ax.set_xlabel('x [nm]', color='#aaa', fontsize=8)
    styled(ax, 'Aerial Image - Hopkins (SVD TCC)')

    fig.suptitle(
        f'Aerial Image Simulation\n'
        f'lam={sys.wavelength:.1f} nm  NA={sys.NA}  sigma={sys.sigma}'
        f'  {sys.N}x{sys.N}  dx={sys.dx:.1f} nm',
        color='white', fontsize=12, y=1.01)

    fname = f'{output_prefix}_aerial_image.png'
    plt.savefig(fname, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    print(f"  Saved: {fname}")
    plt.close()


def plot_linecuts(AI_abbe, AI_hop, sys, output_prefix):
    N  = sys.N
    cx = N // 2
    x  = (np.arange(N) - cx) * sys.dx

    a_cut = AI_abbe[cx, :] / AI_abbe.max()
    h_cut = AI_hop [cx, :] / AI_hop.max()

    fig, ax = plt.subplots(figsize=(9, 4))
    fig.patch.set_facecolor('#0f0f1a')
    ax.set_facecolor('#0f0f1a')
    ax.plot(x, a_cut, color='#ff6b35', lw=2,          label='Abbe')
    ax.plot(x, h_cut, color='#4ecdc4', lw=2, ls='--', label='Hopkins')
    ax.set_xlabel('x [nm]', color='#aaa')
    ax.set_ylabel('Normalised intensity', color='#aaa')
    ax.set_title('Line-cut (centre row)', color='white', fontsize=11)
    ax.legend(facecolor='#1a1a2e', labelcolor='white')
    ax.tick_params(colors='#aaa')
    for sp in ax.spines.values(): sp.set_color('#444')
    ax.grid(color='#333', lw=0.5)

    plt.tight_layout()
    fname = f'{output_prefix}_linecut.png'
    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='#0f0f1a')
    print(f"  Saved: {fname}")
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Aerial image simulation from a source CSV and mask OASIS file.')
    parser.add_argument('--source', required=True,
                        help='Path to source distribution CSV')
    parser.add_argument('--mask',   required=True,
                        help='Path to mask OASIS (.oas) file')
    parser.add_argument('--n_svd',  type=int, default=30,
                        help='Hopkins SVD truncation (default: 30)')
    parser.add_argument('--output', default=None,
                        help='Output filename prefix (default: mask basename)')
    args = parser.parse_args()

    if args.output is None:
        args.output = os.path.splitext(os.path.basename(args.mask))[0]

    print(f"Source : {args.source}")
    print(f"Mask   : {args.mask}")
    print(f"Output : {args.output}_aerial_image.png  /  {args.output}_linecut.png")

    print("\nReading source...")
    source, sys = read_source_csv(args.source)
    print(f"  Non-zero source points : {np.count_nonzero(source)}")
    print(f"  Intensity sum          : {source.sum():.6f}")
    print(f"  lam={sys.wavelength} nm  NA={sys.NA}  "
          f"sigma={sys.sigma}  N={sys.N}  dx={sys.dx} nm")

    print("\nBuilding grids & pupil...")
    X, Y, FX, FY = make_grids(sys)
    pupil = make_pupil(sys, FX, FY, zernike_coeffs={(3, 1): 0.02})

    print("\nReading mask...")
    mask = read_mask_oas(args.mask, sys)
    print(f"  Mask transmission : {mask.mean():.3f}")

    print("\nRunning Abbe simulation...")
    AI_abbe = abbe_simulation(sys, mask, source, pupil, FX, FY)
    print(f"  Intensity range : [{AI_abbe.min():.4f}, {AI_abbe.max():.4f}]")

    print("\nRunning Hopkins simulation...")
    AI_hop = hopkins_simulation(sys, mask, source, pupil, FX, FY, n_svd=args.n_svd)
    print(f"  Intensity range : [{AI_hop.min():.4f}, {AI_hop.max():.4f}]")

    a_n = AI_abbe / AI_abbe.max()
    h_n = AI_hop  / AI_hop.max()
    print(f"\nNorm. RMS(Abbe-Hopkins) = {np.sqrt(np.mean((a_n - h_n)**2)):.5f}")

    print("\nPlotting...")
    plot_results(sys, mask, source, pupil, AI_abbe, AI_hop,
                 X, Y, FX, FY, output_prefix=args.output)
    plot_linecuts(AI_abbe, AI_hop, sys, output_prefix=args.output)

    print("\nDone.")
