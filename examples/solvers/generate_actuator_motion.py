#!/usr/bin/env python3
"""
Generate prescribed circular motion file for the actuator BCP.

Physical setup:
  - Actuator at rest position (8, 0, 0) m
  - Circular motion in the XZ plane: radius R, period T
  - Half-cosine ramp over T_ramp seconds for smooth startup

Motion equations (with ramp):
  theta(t) = omega * t
  ramp(t)  = 0.5 * (1 - cos(pi * t / T_ramp))  for t < T_ramp, else 1.0
  dramp(t) = 0.5 * pi / T_ramp * sin(pi * t / T_ramp)  for t < T_ramp, else 0.0
  d2ramp(t)= 0.5 * (pi/T_ramp)^2 * cos(pi*t/T_ramp) for t < T_ramp, else 0.0

  x(t) = x0 + R * sin(theta) * ramp
  y(t) = 0
  z(t) = z0 + R * (cos(theta) - 1) * ramp

Velocities and accelerations computed via product rule.

Output format (dataActuator_0.dat):
  Line 1: number of time steps
  Lines 2+: t  x y z  vx vy vz  ax ay az
"""
import numpy as np
import os

# ============================================================
# Parameters
# ============================================================
R = 0.5           # Radius of circular motion [m]
T = 10.0          # Period of circular motion [s]
T_ramp = 20.0     # Ramp-up duration [s] (2 full periods)
T_sim = 50.0      # Total simulation time [s]
dt_file = 0.001   # Time resolution for actuator file [s]

x0 = 8.0          # Rest position X [m]
y0 = 0.0          # Rest position Y [m]
z0 = 0.0          # Rest position Z [m]

omega = 2.0 * np.pi / T  # Angular frequency [rad/s]

# Example directories to write the actuator file into
EXAMPLES = ["bdf1", "bdfn_order2", "bdfn_order4", "esdirk46"]

# ============================================================
# Generate time vector
# ============================================================
t = np.arange(0.0, T_sim + 0.5 * dt_file, dt_file)
nt = len(t)

# ============================================================
# Ramp function and its derivatives
# ============================================================
ramp = np.where(t < T_ramp,
                0.5 * (1.0 - np.cos(np.pi * t / T_ramp)),
                1.0)

dramp = np.where(t < T_ramp,
                 0.5 * np.pi / T_ramp * np.sin(np.pi * t / T_ramp),
                 0.0)

d2ramp = np.where(t < T_ramp,
                  0.5 * (np.pi / T_ramp) ** 2 * np.cos(np.pi * t / T_ramp),
                  0.0)

# ============================================================
# Angular quantities
# ============================================================
theta = omega * t
sin_th = np.sin(theta)
cos_th = np.cos(theta)

# ============================================================
# Position
# ============================================================
x = x0 + R * sin_th * ramp
y = np.full_like(t, y0)
z = z0 + R * (cos_th - 1.0) * ramp

# ============================================================
# Velocity (product rule: d/dt [f(t)*ramp(t)])
# ============================================================
# dx/dt = R * [omega*cos(theta)*ramp + sin(theta)*dramp]
vx = R * (omega * cos_th * ramp + sin_th * dramp)
vy = np.zeros_like(t)
# dz/dt = R * [-omega*sin(theta)*ramp + (cos(theta)-1)*dramp]
vz = R * (-omega * sin_th * ramp + (cos_th - 1.0) * dramp)

# ============================================================
# Acceleration (product rule on velocity terms)
# ============================================================
# d2x/dt2 = R * [-omega^2*sin(theta)*ramp + 2*omega*cos(theta)*dramp + sin(theta)*d2ramp]
ax = R * (-omega**2 * sin_th * ramp + 2.0 * omega * cos_th * dramp + sin_th * d2ramp)
ay = np.zeros_like(t)
# d2z/dt2 = R * [-omega^2*cos(theta)*ramp - 2*omega*sin(theta)*dramp + (cos(theta)-1)*d2ramp]
az = R * (-omega**2 * cos_th * ramp - 2.0 * omega * sin_th * dramp + (cos_th - 1.0) * d2ramp)

# ============================================================
# Write output files
# ============================================================
script_dir = os.path.dirname(os.path.abspath(__file__))

for example in EXAMPLES:
    output_dir = os.path.join(script_dir, example, "input")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "dataActuator_0.dat")

    with open(output_file, "w") as f:
        f.write(f"{nt}\n")
        for i in range(nt):
            f.write(f"{t[i]:.4f} "
                    f"{x[i]:.10f} {y[i]:.10f} {z[i]:.10f} "
                    f"{vx[i]:.10f} {vy[i]:.10f} {vz[i]:.10f} "
                    f"{ax[i]:.10f} {ay[i]:.10f} {az[i]:.10f}\n")

    print(f"Written {output_file} ({nt} time steps)")

print("Done.")
