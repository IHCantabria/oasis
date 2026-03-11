"""
Plot results for the large_rotation_example.
Compares free-decay motions between simplified (rotSimpFlag=1) and
full rotation dynamics (rotSimpFlag=0).
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Configuration
output_simplified = os.path.join("case_simplified", "output")
output_full = os.path.join("case_full", "output")
body_id = 0
dof_labels = ["Surge (x)", "Sway (y)", "Heave (z)", "Roll", "Pitch", "Yaw"]
dof_units = ["m", "m", "m", "rad", "rad", "rad"]

fig, axes = plt.subplots(3, 2, figsize=(14, 10), sharex=True)
fig.suptitle("Large Rotation Example - Simplified vs Full Rotation Dynamics\n(Free decay: heave + pitch initial displacement)", fontsize=13)

for i in range(6):
    row = i % 3
    col = i // 3
    ax = axes[row, col]

    # Simplified
    f_simp = os.path.join(output_simplified, f"DOF_{i+1}_Body_{body_id}.txt")
    if os.path.exists(f_simp):
        data = np.loadtxt(f_simp)
        ax.plot(data[:, 0], data[:, 1], "b-", linewidth=1.0, label="Simplified (rotSimpFlag=1)")

    # Full
    f_full = os.path.join(output_full, f"DOF_{i+1}_Body_{body_id}.txt")
    if os.path.exists(f_full):
        data = np.loadtxt(f_full)
        ax.plot(data[:, 0], data[:, 1], "r--", linewidth=1.0, label="Full (rotSimpFlag=0)")

    ax.set_ylabel(f"{dof_labels[i]} [{dof_units[i]}]")
    ax.set_title(f"DOF {i+1}: {dof_labels[i]}")
    ax.grid(True, alpha=0.3)
    if i == 0:
        ax.legend(fontsize=8)

for ax in axes[-1, :]:
    ax.set_xlabel("Time [s]")

plt.tight_layout()
plt.savefig("large_rotation_results.png", dpi=150)
print("Saved large_rotation_results.png")
plt.show()
