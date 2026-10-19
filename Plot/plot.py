#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt


if len(sys.argv) != 2:
    print(f"Usage: python {sys.argv[0]} <filename>")
    sys.exit(1)

filename = sys.argv[1]

# Read the whole file
with open(filename, "r") as f:
    lines = f.readlines()

# First line: dimensions
nx, ny = map(int, lines[0].split())

# The file contains:
#   nx+1 lines -> first scalar field
#   "vector"
#   nx+1 lines -> second field
#
# Find "vector"
vector_line = next(i for i, line in enumerate(lines) if line.strip() == "vector")

# Read first field
field1 = np.array([
    [float(x) for x in lines[i].split()]
    for i in range(1, vector_line)
])

# Read second field
field2 = np.array([
    [float(x) for x in lines[i].split()]
    for i in range(vector_line + 1, len(lines))
])

print("Field 1 shape:", field1.shape)
print("Field 2 shape:", field2.shape)

# Coordinates
x = np.arange(field1.shape[1])
y = np.arange(field1.shape[0])

# Plot
fig, ax = plt.subplots()

im = ax.pcolormesh(x, y, field1, shading="auto")
fig.colorbar(im, ax=ax, label="Field 1")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")

plt.tight_layout()
plt.show()
