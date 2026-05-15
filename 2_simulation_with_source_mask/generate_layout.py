"""
Generate multi_L-shaped.oas — 8 L-shaped polygons with random orientations.

L-shape geometry (before rotation), centred at origin:
    arm width  w = 32 nm
    arm length h = 96 nm  (3 × w)

    (-48, 48) ──── (-16, 48)
        |               |
        |               |
    (-48,-16) ── (-16,-16) ── (48,-16)
        |                         |
    (-48,-48) ──────────────── (48,-48)

Each of the 8 instances is rotated by a random angle and placed on a
4 × 2 grid (pitch = 220 nm) so polygons do not overlap.

Usage:
    python generate_layout.py [--seed SEED] [--output PATH]
"""

import argparse
import math
import random
import numpy as np

def l_shape_vertices(w=32.0, h=96.0):
    """Return the 6 vertices of an L-shape centred at the origin."""
    hw = h / 2
    return np.array([
        [-hw,    -hw   ],
        [ hw,    -hw   ],
        [ hw,    -hw+w ],
        [-hw+w,  -hw+w ],
        [-hw+w,   hw   ],
        [-hw,     hw   ],
    ], dtype=float)


def rotate(pts, angle_deg):
    """Rotate 2-D points around the origin by angle_deg (degrees)."""
    a = math.radians(angle_deg)
    c, s = math.cos(a), math.sin(a)
    R = np.array([[c, -s], [s, c]])
    return pts @ R.T


def translate(pts, dx, dy):
    return pts + np.array([dx, dy])


def write_oasis(polygons, filepath, layer=1, datatype=0):
    try:
        import gdstk
    except ImportError:
        raise ImportError("gdstk is required: pip install gdstk")

    lib  = gdstk.Library()
    cell = lib.new_cell("TOP")

    for pts in polygons:
        poly = gdstk.Polygon(pts, layer=layer, datatype=datatype)
        cell.add(poly)

    lib.write_oas(filepath)
    import os
    print(f"Saved: {filepath}  ({os.path.getsize(filepath)/1024:.1f} kB, "
          f"{len(polygons)} polygons)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create multi_L-shaped.oas with 8 randomly oriented L polygons.")
    parser.add_argument("--seed",   type=int,   default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("--output", type=str,   default="multi_L-shaped.oas",
                        help="Output OASIS file path")
    args = parser.parse_args()

    random.seed(args.seed)

    # 4 columns × 2 rows, 220 nm pitch.
    # Offset by 150 nm so rotated L-shapes (max radius ~68 nm) stay inside
    # the simulation grid [0, 1024] nm with comfortable margin.
    cols, rows, pitch, offset = 4, 2, 220.0, 150.0
    positions = [
        (offset + c * pitch, offset + r * pitch)
        for r in range(rows)
        for c in range(cols)
    ]  # 8 positions total

    angles = [random.uniform(0, 360) for _ in range(8)]

    print("L-shape orientations:")
    polygons = []
    for i, ((cx, cy), angle) in enumerate(zip(positions, angles)):
        pts = l_shape_vertices(w=32.0, h=96.0)
        pts = rotate(pts, angle)
        pts = translate(pts, cx, cy)
        polygons.append(pts)
        print(f"  [{i+1}] centre=({cx:.0f}, {cy:.0f}) nm   rotation={angle:.1f}°")

    write_oasis(polygons, args.output)
