#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt


def read_mco(filename):
    filename = './Output/Output_2D/' + filename + '.mco'

    with open(filename, "r") as f:
        # First line: 49 49
        nx, ny = map(int, f.readline().split())

        # Read the remaining values
        data = np.loadtxt(f)

    data = data.reshape((ny + 1, nx + 1))

    return data


# --------------------------------------------------
# Command line arguments
# --------------------------------------------------

if len(sys.argv) != 2:
    print(f"Usage: python {sys.argv[0]} <Qty>")
    sys.exit(1)

Qty = sys.argv[1]

# The three coordinate planes to plot. Component files are assumed to be
# named "<Qty><component>", e.g. "Bx", "By", "Bz" for Qty="B" — matching
# the single-pair script's convention where argv gave the two component
# filenames directly.
planes = [("x", "y"), ("x", "z"), ("y", "z")]


# --------------------------------------------------
# Figure with one subplot per plane
# --------------------------------------------------

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for ax, (direction1, direction2) in zip(axes, planes):

    # ----------------------------------------------
    # Read fields for this plane's two components
    # ----------------------------------------------
    comp1 = read_mco(Qty + direction1 + '_' + direction1 + direction2)
    comp2 = read_mco(Qty + direction2 + '_' + direction1 + direction2)

    if comp1.shape != comp2.shape:
        raise ValueError(
            f"Fields have different shapes: {comp1.shape} vs {comp2.shape}"
        )

    # ----------------------------------------------
    # Grid
    # ----------------------------------------------
    ny, nx = comp1.shape

    y = np.arange(ny)
    z = np.arange(nx)

    Y, Z = np.meshgrid(y, z, indexing="ij")
    
    x_idx = np.arange(nx)
    y_idx = np.arange(ny)
    Xg, Yg = np.meshgrid(x_idx, y_idx) 

    # ----------------------------------------------
    # Vector magnitude
    # ----------------------------------------------
    mag = np.sqrt(comp1**2 + comp2**2)

    # ----------------------------------------------
    # Plot: magnitude as background, plain arrows on top
    # ----------------------------------------------


    mesh = ax.pcolormesh(
        Xg, Yg,
        mag,
        cmap="rainbow",
        alpha=0.45,
        shading="auto",
    )

    cbar = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$|\mathbf{B}|$")


    step = 4  # Plot every 4th grid point

    ax.quiver(
        Xg[1::step, 1::step],
        Yg[1::step, 1::step],
        comp1[1::step, 1::step],
        comp2[1::step, 1::step],
        color="black",
        angles="xy",
        scale_units="xy",
        scale=0.2*np.max(mag),  # Smaller value = longer arrows
        width=0.004,
    )

    ax.set_xlabel(direction1)
    ax.set_ylabel(direction2)
    ax.set_title(f"{Qty}_{direction1}{direction2}")
    ax.set_aspect("equal")

   

plt.tight_layout()
plt.show()