#!/usr/bin/env python3
"""
Plot an ASCII .mco plane file written by mod_io_legacy_mco.f90.

Usage:
    python plot_mco.py path/to/file.mco
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def read_mco_ascii(path: Path) -> np.ndarray:
    """
    Read ASCII mco: first line 'n1 n2', then n2 rows with n1 values.
    Returns array with shape (n2, n1).
    """
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()

        if len(header) != 2:
            raise ValueError(
                f"Bad header in {path.name}: expected 'n1 n2', got {header}"
            )

        n1 = int(header[0]) + 1
        n2 = int(header[1]) + 1

        data = []
        for _ in range(n2):
            row = f.readline().strip().split()
            if len(row) != n1:
                raise ValueError(
                    f"Row has {len(row)} values but expected {n1}"
                )
            data.append([float(x) for x in row])

    return np.array(data, dtype=float)


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: python plot_mco.py <file.mco>")
        return 1

    path = Path(sys.argv[1])

    if not path.is_file():
        print(f"Error: file not found: {path}")
        return 1

    try:
        arr = read_mco_ascii(path)
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        return 1

    # --- Plot ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # 2D plot
    im = ax1.imshow(arr, origin="lower", aspect="auto")
    plt.colorbar(im, ax=ax1, label=r"$\phi$ (V)")
    ax1.set_xlabel("x cell")
    ax1.set_ylabel("z cell")

    # 1D slice
    #ax2.plot(arr[16, :], label="LePIC potential", linewidth=3)
    ax2.plot()
    plt.tight_layout()
    outname = path.with_suffix(".png")
    plt.savefig(outname, dpi=200, bbox_inches="tight")
    print(f"Saved figure to {outname}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
