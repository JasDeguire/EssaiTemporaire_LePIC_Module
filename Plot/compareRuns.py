#!/usr/bin/env python3
"""
Compare Output/Output_2D/itXXXX_name_dim.mco between two iterations.

Usage:
    python compareRuns.py name it1 it2

Example:
    python compareRuns.py phi 1001 1002
        -> compares it1001_phi_{xy,xz,yz}.mco against it1002_phi_{xy,xz,yz}.mco

Both files are read from Output/Output_2D under the current directory.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def read_mco_ascii(path: Path) -> np.ndarray:
    data = []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()

        if len(header) != 2:
            raise ValueError(f"Bad header in {path}: expected 'n1 n2', got {header}")

        for line in f:
            row = line.strip().split()
            if row:
                data.append([float(x) for x in row])

    if not data:
        raise ValueError(f"No data found in {path}")

    ncols = len(data[0])

    for i, row in enumerate(data):
        if len(row) != ncols:
            raise ValueError(
                f"{path}: row {i + 1} has {len(row)} values, expected {ncols}"
            )

    return np.array(data, dtype=float)


def resample_to_shape(arr: np.ndarray, target_shape: tuple[int, int]) -> np.ndarray:
    ny_old, nx_old = arr.shape
    ny_new, nx_new = target_shape

    if (ny_old, nx_old) == (ny_new, nx_new):
        return arr

    x_old = np.linspace(0.0, 1.0, nx_old)
    y_old = np.linspace(0.0, 1.0, ny_old)

    x_new = np.linspace(0.0, 1.0, nx_new)
    y_new = np.linspace(0.0, 1.0, ny_new)

    tmp = np.empty((ny_old, nx_new))

    for j in range(ny_old):
        tmp[j, :] = np.interp(x_new, x_old, arr[j, :])

    out = np.empty((ny_new, nx_new))

    for i in range(nx_new):
        out[:, i] = np.interp(y_new, y_old, tmp[:, i])

    return out


def axis_labels(dim: str) -> tuple[str, str]:
    if dim == "xy":
        return "normalized x", "normalized y"
    if dim == "xz":
        return "normalized x", "normalized z"
    if dim == "yz":
        return "normalized y", "normalized z"

    return "normalized axis 1", "normalized axis 2"


def compare_one_dim(
    name: str,
    dim: str,
    it1: int,
    it2: int,
    root_dir: Path,
) -> None:

    path1 = root_dir / "Output" / "Output_2D" / f"it{it1}_{name}_{dim}.mco"
    path2 = root_dir / "Output" / "Output_2D" / f"it{it2}_{name}_{dim}.mco"

    if not path1.exists():
        raise FileNotFoundError(f"Missing: {path1}")

    if not path2.exists():
        raise FileNotFoundError(f"Missing: {path2}")

    arr1 = read_mco_ascii(path1)
    arr2 = read_mco_ascii(path2)

    print(f"\n{dim}:")
    print(f"  {path1.name} shape = {arr1.shape}")
    print(f"  {path2.name} shape = {arr2.shape}")

    arr2_interp = resample_to_shape(arr2, arr1.shape)
    diff = arr2_interp - arr1

    print(f"  comparison shape = {diff.shape}")
    print(f"  max abs diff     = {np.max(np.abs(diff)):.6e}")
    print(f"  mean diff        = {np.mean(diff):.6e}")

    xlabel, ylabel = axis_labels(dim)

    fig, axes = plt.subplots(1, 4, figsize=(24, 5))

    extent = (0.0, 1.0, 0.0, 1.0)
    maximum = max(np.max(arr1), np.max(arr2_interp))
    minimum = min(np.min(arr1), np.min(arr2_interp))
    im1 = axes[0].imshow(arr1, origin="lower", aspect="auto", extent=extent, vmin=minimum, vmax=maximum, cmap='Spectral')
    plt.colorbar(im1, ax=axes[0])
    axes[0].set_title(f"it{it1}\n{path1.name}")
    axes[0].set_xlabel(xlabel)
    axes[0].set_ylabel(ylabel)

    im2 = axes[1].imshow(arr2_interp, origin="lower", aspect="auto", extent=extent, vmin=minimum, vmax=maximum, cmap='Spectral')
    plt.colorbar(im2, ax=axes[1])
    axes[1].set_title(f"it{it2}\n{path2.name}")
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel(ylabel)

    im3 = axes[2].imshow(diff, origin="lower", aspect="auto", extent=extent)
    plt.colorbar(im3, ax=axes[2])
    axes[2].set_title(f"it{it2} - it{it1}")
    axes[2].set_xlabel(xlabel)
    axes[2].set_ylabel(ylabel)

    avg1 = np.mean(arr1, axis=0)
    avg2 = np.mean(arr2_interp, axis=0)
    x = np.linspace(0.0, 1.0, avg1.size)

    axes[3].plot(x, avg1, label=f"it{it1}", linewidth=1)
    axes[3].plot(x, avg2, label=f"it{it2}", linewidth=1)

    axes[3].set_title(f"Average profiles, {dim}")
    axes[3].set_xlabel(xlabel)
    axes[3].set_ylabel(name)
    axes[3].legend()
    axes[3].grid(True, alpha=0.3)

    fig.suptitle(f"{name}_{dim} — it{it1} vs it{it2}", fontsize=14)

    plt.tight_layout()

    output_file = root_dir / f"compareRuns_{name}_{dim}_it{it1}_vs_it{it2}.png"
    plt.savefig(output_file, dpi=200)
    plt.close(fig)

    print(f"  saved: {output_file}")


def main() -> int:
    if len(sys.argv) != 4:
        print("Usage:")
        print("    python compareRuns.py name it1 it2")
        print("")
        print("Example:")
        print("    python compareRuns.py phi 1001 1002")
        print("        -> compares it1001_phi_*.mco against it1002_phi_*.mco")
        return 1

    name = sys.argv[1]

    try:
        it1 = int(sys.argv[2])
        it2 = int(sys.argv[3])
    except ValueError:
        print(f"it1/it2 must be integers, got '{sys.argv[2]}' and '{sys.argv[3]}'")
        return 1

    root_dir = Path.cwd()

    for dim in ("xy", "xz", "yz"):
        try:
            compare_one_dim(name, dim, it1, it2, root_dir)
        except Exception as e:
            print(f"Error for dim={dim}: {e}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
