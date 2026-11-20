#!/usr/bin/env python3
"""
Compare legacy DATA/DATA_2D/name_dim.mco
with modular Output/Output_2D/itXXXX_name_dim.mco.

Handles different shapes by using normalized coordinates.
The legacy array is interpolated onto the modular grid before subtraction.

Usage:
    python compare_mco.py name
    python compare_mco.py name iteration
    python compare_mco.py name start:stop:step

Examples:
    python compare_mco.py phi
    python compare_mco.py phi 2001
    python compare_mco.py phi 120001:150001:1000

Average mode only averages modular files.
Legacy file is still read once:
    DATA/DATA_2D/name_dim.mco
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


def parse_average_range(spec: str) -> list[int]:
    try:
        start_s, stop_s, step_s = spec.split(":")
        start = int(start_s)
        stop = int(stop_s)
        step = int(step_s)
    except ValueError:
        raise ValueError(
            "Average range must have format start:stop:step, "
            "for example 120001:150001:1000"
        )

    if step <= 0:
        raise ValueError("Step must be positive")

    return list(range(start, stop + 1, step))


def apply_name_scaling(name: str, arr: np.ndarray) -> np.ndarray:
    arr = arr.copy()



    return arr


def read_modular_average(
    name: str,
    dim: str,
    iterations: list[int],
    root_dir: Path,
) -> tuple[np.ndarray, list[Path]]:

    arrays = []
    files = []

    for it in iterations:
        path = root_dir / "Output" / "Output_2D" / f"it{it}_{name}_{dim}.mco"

        if not path.exists():
            print(f"Missing: {path.name}")
            continue

        try:
            arr = read_mco_ascii(path)
            arrays.append(arr)
            files.append(path)
            print(f"Read OK: {path.name}")
        except Exception as e:
            print(f"Skipping {path.name}: {e}")

    if not arrays:
        raise RuntimeError(f"No valid modular files found for {name}_{dim}")

    shape0 = arrays[0].shape

    for path, arr in zip(files, arrays):
        if arr.shape != shape0:
            raise ValueError(
                f"Shape mismatch in {path.name}: got {arr.shape}, expected {shape0}"
            )

    avg = np.mean(arrays, axis=0)

    return avg, files


def compare_arrays(
    name: str,
    dim: str,
    arr_legacy: np.ndarray,
    arr_modular_raw: np.ndarray,
    root_dir: Path,
    modular_label: str,
    output_label: str,
) -> None:

    print(f"\n{dim}:")
    print(f"  legacy shape  = {arr_legacy.shape}")
    print(f"  modular shape = {arr_modular_raw.shape}")

    arr_modular = apply_name_scaling(name, arr_modular_raw)

    arr_legacy_interp = resample_to_shape(arr_legacy, arr_modular.shape)
    diff = arr_modular - arr_legacy_interp

    print(f"  comparison shape = {diff.shape}")
    print(f"  max abs diff     = {np.max(np.abs(diff)):.6e}")
    print(f"  mean diff        = {np.mean(diff):.6e}")

    xlabel, ylabel = axis_labels(dim)

    fig, axes = plt.subplots(1, 4, figsize=(24, 5))

    extent = (0.0, 1.0, 0.0, 1.0)
    maximum = max(np.max(arr_legacy), np.max(arr_modular))
    minimum = min(np.min(arr_legacy), np.min(arr_modular))
    im1 = axes[0].imshow(
        arr_legacy,
        origin="lower",
        aspect="auto",
        extent=extent,
        vmin=minimum,
        vmax=maximum)

    plt.colorbar(im1, ax=axes[0])
    axes[0].set_title(f"Legacy\n{name}_{dim}.mco\nshape={arr_legacy.shape}")
    axes[0].set_xlabel(xlabel)
    axes[0].set_ylabel(ylabel)

    im2 = axes[1].imshow(
        arr_modular,
        origin="lower",
        aspect="auto",
        extent=extent,
        vmin=minimum,
        vmax=maximum)
    plt.colorbar(im2, ax=axes[1])
    axes[1].set_title(f"Modular\n{modular_label}\nshape={arr_modular.shape}")
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel(ylabel)

    im3 = axes[2].imshow(
        diff,
        origin="lower",
        aspect="auto",
        extent=extent,
    )
    plt.colorbar(im3, ax=axes[2])
    axes[2].set_title("Modular - legacy interpolated")
    axes[2].set_xlabel(xlabel)
    axes[2].set_ylabel(ylabel)

    avg_legacy = np.mean(arr_legacy, axis=0)
    avg_modular = np.mean(arr_modular, axis=0)

    x_legacy = np.linspace(0.0, 1.0, avg_legacy.size)
    x_modular = np.linspace(0.0, 1.0, avg_modular.size)

    avg_legacy_interp = np.interp(x_modular, x_legacy, avg_legacy)

    axes[3].plot(
        x_modular,
        avg_legacy_interp,
        label="legacy interpolated",
        linewidth=3,
    )
    axes[3].plot(
        x_modular,
        avg_modular,
        label="modular",
        linewidth=3,
    )

    axes[3].set_title(f"Average profiles, {dim}")
    axes[3].set_xlabel(xlabel)
    axes[3].set_ylabel(name)
    axes[3].legend()
    axes[3].grid(True, alpha=0.3)

    fig.suptitle(f"{name}_{dim} — normalized comparison", fontsize=14)

    plt.tight_layout()

    output_file = root_dir / f"compare_{output_label}_{name}_{dim}.png"
    plt.savefig(output_file, dpi=200)
    plt.close(fig)

    print(f"  saved: {output_file}")


def compare_one_dim(
    name: str,
    dim: str,
    iteration: int,
    root_dir: Path,
) -> None:

    legacy_path = root_dir / "DATA" / "DATA_2D" / f"{name}_{dim}.mco"
    modular_path = root_dir / "Output" / "Output_2D" / f"it{iteration}_{name}_{dim}.mco"

    if not legacy_path.exists():
        raise FileNotFoundError(f"Missing legacy file: {legacy_path}")

    if not modular_path.exists():
        raise FileNotFoundError(f"Missing modular file: {modular_path}")

    arr_legacy = read_mco_ascii(legacy_path)
    arr_modular = read_mco_ascii(modular_path)

    compare_arrays(
        name=name,
        dim=dim,
        arr_legacy=arr_legacy,
        arr_modular_raw=arr_modular,
        root_dir=root_dir,
        modular_label=modular_path.name,
        output_label=f"it{iteration}",
    )


def compare_average_dim(
    name: str,
    dim: str,
    iterations: list[int],
    root_dir: Path,
) -> None:

    legacy_path = root_dir / "DATA" / "DATA_2D" / f"{name}_{dim}.mco"

    if not legacy_path.exists():
        raise FileNotFoundError(f"Missing legacy file: {legacy_path}")

    arr_legacy = read_mco_ascii(legacy_path)

    arr_modular_avg, files = read_modular_average(
        name=name,
        dim=dim,
        iterations=iterations,
        root_dir=root_dir,
    )

    print(f"\nAveraged {len(files)} modular files for {dim}")

    modular_label = f"average it{iterations[0]} to it{iterations[-1]}"
    output_label = f"avg_it{iterations[0]}_to_it{iterations[-1]}"

    compare_arrays(
        name=name,
        dim=dim,
        arr_legacy=arr_legacy,
        arr_modular_raw=arr_modular_avg,
        root_dir=root_dir,
        modular_label=modular_label,
        output_label=output_label,
    )


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage:")
        print("    python compare_mco.py name")
        print("    python compare_mco.py name iteration")
        print("    python compare_mco.py name start:stop:step")
        print("")
        print("Examples:")
        print("    python compare_mco.py phi")
        print("    python compare_mco.py phi 2001")
        print("    python compare_mco.py phi 120001:150001:1000")
        return 1

    name = sys.argv[1]
    root_dir = Path.cwd()

    if len(sys.argv) == 2:
        iteration = 110001

        for dim in ("xy", "xz", "yz"):
            try:
                compare_one_dim(name, dim, iteration, root_dir)
            except Exception as e:
                print(f"Error for dim={dim}: {e}", file=sys.stderr)

        return 0

    arg = sys.argv[2]

    if ":" in arg:
        try:
            iterations = parse_average_range(arg)
        except Exception as e:
            print(f"Bad averaging range: {e}")
            return 1

        print(f"Averaging modular iterations from {iterations[0]} to {iterations[-1]}")

        for dim in ("xy", "xz", "yz"):
            try:
                compare_average_dim(name, dim, iterations, root_dir)
            except Exception as e:
                print(f"Error for dim={dim}: {e}", file=sys.stderr)

    else:
        try:
            iteration = int(arg)
        except ValueError:
            print(
                "Second argument must be either:\n"
                "  iteration number\n"
                "or\n"
                "  start:stop:step"
            )
            return 1

        for dim in ("xy", "xz", "yz"):
            try:
                compare_one_dim(name, dim, iteration, root_dir)
            except Exception as e:
                print(f"Error for dim={dim}: {e}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())