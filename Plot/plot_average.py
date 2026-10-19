#!/usr/bin/env python3
"""
Plot modular Output/Output_2D/itXXXX_name_dim.mco files.

Usage:
    python plot_mco_modular.py name
    python plot_mco_modular.py name iteration
    python plot_mco_modular.py name start:stop:step

Examples:
    python plot_mco_modular.py phi
    python plot_mco_modular.py phi 2001
    python plot_mco_modular.py n1 120001:150001:1000

The averaging mode:
    start:stop:step

Example:
    120001:150001:1000

will average:
    it120001_name_dim.mco
    it121001_name_dim.mco
    ...
    it150001_name_dim.mco
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


# ---------------------------------
# READ MCO
# ---------------------------------
def read_mco_ascii(path: Path) -> np.ndarray:
    data = []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()

        if len(header) != 2:
            raise ValueError(
                f"Bad header in {path}: expected 'n1 n2', got {header}"
            )

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
                f"{path}: row {i + 1} has {len(row)} values, "
                f"expected {ncols}"
            )

    return np.array(data, dtype=float)


# ---------------------------------
# WRITE MCO
# ---------------------------------
def write_mco_ascii(path: Path, arr: np.ndarray) -> None:
    n2, n1 = arr.shape

    with path.open("w", encoding="utf-8") as f:
        f.write(f"{n1-1} {n2-1}\n")

        for j in range(n2):
            f.write(
                " ".join(f"{v:.6e}" for v in arr[j, :]) + "\n"
            )


# ---------------------------------
# AXIS LABELS
# ---------------------------------
def axis_labels(dim: str) -> tuple[str, str]:

    if dim == "xy":
        return "x index", "y index"

    if dim == "xz":
        return "x index", "z index"

    if dim == "yz":
        return "y index", "z index"

    return "axis 1 index", "axis 2 index"


# ---------------------------------
# AVERAGE RANGE
# ---------------------------------
def parse_average_range(spec: str) -> list[int]:

    try:
        start_s, stop_s, step_s = spec.split(":")

        start = int(start_s)
        stop = int(stop_s)
        step = int(step_s)

    except ValueError:
        raise ValueError(
            "Average range must have format:\n"
            "start:stop:step\n"
            "Example:\n"
            "120001:150001:1000"
        )

    if step <= 0:
        raise ValueError("Step must be positive")

    return list(range(start, stop + 1, step))


# ---------------------------------
# SPECIAL SCALING
# ---------------------------------
def apply_name_scaling(
    name: str,
    arr: np.ndarray
) -> np.ndarray:

    arr = arr.copy()



    return arr


# ---------------------------------
# READ AVERAGE
# ---------------------------------
def read_average(
    name: str,
    dim: str,
    iterations: list[int],
    root_dir: Path,
) -> tuple[np.ndarray, list[Path]]:

    files = []
    arrays = []

    for it in iterations:

        path = (
            root_dir
            / "Output"
            / "Output_2D"
            / f"it{it}_{name}_{dim}.mco"
        )

        if not path.exists():
            print(f"Missing: {path.name}")
            continue

        try:
            arr = read_mco_ascii(path)

            files.append(path)
            arrays.append(arr)

            print(f"Read OK: {path.name}")

        except Exception as e:
            print(f"Skipping {path.name}: {e}")

    if not arrays:
        raise RuntimeError(
            f"No valid files found for {name}_{dim}"
        )

    shape0 = arrays[0].shape

    for arr in arrays:
        if arr.shape != shape0:
            raise ValueError(
                "Shape mismatch between files."
            )

    avg = np.mean(arrays, axis=0)

    return avg, files


# ---------------------------------
# PLOT ARRAY
# ---------------------------------
def plot_array(
    arr: np.ndarray,
    name: str,
    dim: str,
    root_dir: Path,
    title_prefix: str,
    output_prefix: str,
) -> None:

    arr = apply_name_scaling(name, arr)

    print(f"\n{dim}:")
    print(f"  shape = {arr.shape}")
    print(f"  min   = {np.min(arr):.6e}")
    print(f"  max   = {np.max(arr):.6e}")
    print(f"  mean  = {np.mean(arr):.6e}")

    xlabel, ylabel = axis_labels(dim)

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(14, 5)
    )

    # -----------------------------
    # 2D MAP
    # -----------------------------
    if(name[0] == 'n'):
        np.where((arr < np.max(arr)/10**3) & (arr > np.max(arr)/10**6) , arr, np.max(arr)/10**3)
        im = axes[0].imshow(
            arr,
            origin="lower",
            cmap="Spectral_r",
            aspect="auto",
            norm=LogNorm(vmin=max(np.max(arr)/10**3,np.min(arr)), vmax=np.max(arr))
        )
        axes[1].set_yscale('log')

    else:
        im = axes[0].imshow(
            arr,
            origin="lower",
            cmap="Spectral_r",
            aspect="auto",
        )

    plt.colorbar(im, ax=axes[0])

    axes[0].set_title(
        f"{title_prefix}_{name}_{dim}\n"
        f"shape={arr.shape}"
        )

    axes[0].set_xlabel(xlabel)
    axes[0].set_ylabel(ylabel)

    # -----------------------------
    # 1D PROFILE
    # -----------------------------
    avg_profile = np.mean(arr, axis=0)

    x = np.arange(avg_profile.size)

    axes[1].plot(
        x,
        avg_profile,
        linewidth=3,
    )

    axes[1].set_title(
        f"Average profile, {dim}"
    )

    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel(name)

    axes[1].grid(
        True,
        alpha=0.3
    )

    fig.suptitle(
        f"{title_prefix} {name}_{dim}",
        fontsize=14
    )

    plt.tight_layout()

    output_file = (
        root_dir
        / f"{output_prefix}_{name}_{dim}.png"
    )

    plt.savefig(
        output_file,
        dpi=200
    )

    plt.close(fig)

    print(f"Saved plot: {output_file}")


# ---------------------------------
# SINGLE ITERATION
# ---------------------------------
def plot_one_dim(
    name: str,
    dim: str,
    iteration: int,
    root_dir: Path,
) -> None:

    path = (
        root_dir
        / "Output"
        / "Output_2D"
        / f"it{iteration}_{name}_{dim}.mco"
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Missing modular file: {path}"
        )

    arr = read_mco_ascii(path)

    print(f"\nUsing file:")
    print(f"  {path.name}")

    plot_array(
        arr=arr,
        name=name,
        dim=dim,
        root_dir=root_dir,
        title_prefix=f"it{iteration}",
        output_prefix=f"plot_it{iteration}",
    )


# ---------------------------------
# AVERAGE MODE
# ---------------------------------
def plot_average_dim(
    name: str,
    dim: str,
    iterations: list[int],
    root_dir: Path,
) -> None:

    avg, files = read_average(
        name,
        dim,
        iterations,
        root_dir,
    )

    print(f"\nAveraged {len(files)} files.")

    avg_file = (
        root_dir
        / f"avg_{name}_{dim}.mco"
    )

    write_mco_ascii(
        avg_file,
        avg
    )

    print(f"Saved average file:")
    print(f"  {avg_file}")

    plot_array(
        arr=avg,
        name=name,
        dim=dim,
        root_dir=root_dir,
        title_prefix="average",
        output_prefix="plot_avg",
    )


# ---------------------------------
# MAIN
# ---------------------------------
def main() -> int:

    if len(sys.argv) < 2:

        print("Usage:")
        print("    python plot_mco_modular.py name")
        print("    python plot_mco_modular.py name iteration")
        print("    python plot_mco_modular.py name start:stop:step")
        print("")
        print("Examples:")
        print("    python plot_mco_modular.py phi")
        print("    python plot_mco_modular.py phi 2001")
        print(
            "    python plot_mco_modular.py "
            "n1 120001:150001:1000"
        )

        return 1

    name = sys.argv[1]

    root_dir = Path.cwd()

    # ---------------------------------
    # DEFAULT MODE
    # ---------------------------------
    if len(sys.argv) == 2:

        iteration = 118001

        for dim in ("xy", "xz", "yz"):

            try:
                plot_one_dim(
                    name,
                    dim,
                    iteration,
                    root_dir
                )

            except Exception as e:
                print(
                    f"Error for dim={dim}: {e}",
                    file=sys.stderr
                )

        return 0

    # ---------------------------------
    # SECOND ARGUMENT
    # ---------------------------------
    arg = sys.argv[2]

    # ---------------------------------
    # AVERAGE MODE
    # ---------------------------------
    if ":" in arg:

        try:
            iterations = parse_average_range(arg)

        except Exception as e:
            print(f"Bad averaging range: {e}")
            return 1

        print(
            f"\nAveraging iterations from "
            f"{iterations[0]} "
            f"to {iterations[-1]}"
        )

        for dim in ("xy", "xz", "yz"):

            try:
                plot_average_dim(
                    name,
                    dim,
                    iterations,
                    root_dir
                )

            except Exception as e:
                print(
                    f"Error for dim={dim}: {e}",
                    file=sys.stderr
                )

    # ---------------------------------
    # SINGLE ITERATION MODE
    # ---------------------------------
    else:

        try:
            iteration = int(arg)

        except ValueError:

            print(
                "Second argument must be:\n"
                "  iteration number\n"
                "or\n"
                "  start:stop:step"
            )

            return 1

        for dim in ("xy", "xz", "yz"):

            try:
                plot_one_dim(
                    name,
                    dim,
                    iteration,
                    root_dir
                )

            except Exception as e:
                print(
                    f"Error for dim={dim}: {e}",
                    file=sys.stderr
                )

    return 0


# ---------------------------------
# ENTRY POINT
# ---------------------------------
if __name__ == "__main__":
    raise SystemExit(main())