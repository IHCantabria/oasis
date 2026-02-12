# -*- coding: utf-8 -*-

import numpy as np


def m1(t):

    r = 0.2
    T = 1.25
    w = -2.0 * np.pi / T
    if t > 5:
        x = 0.0 + r * np.cos(w * (t - 5.0) - 0.5 * np.pi)
        z = 0.0 + r * np.sin(w * (t - 5.0) - 0.5 * np.pi)
        vx = -r * w * np.sin(w * (t - 5.0) - 0.5 * np.pi)
        vz = r * w * np.cos(w * (t - 5.0) - 0.5 * np.pi)
        ax = -r * w * w * np.cos(w * (t - 5.0) - 0.5 * np.pi)
        az = -r * w * w * np.sin(w * (t - 5.0) - 0.5 * np.pi)
    else:
        x = 0.0
        z = 0.0
        vx = 0.0
        vz = 0.0
        ax = 0.0
        az = 0.0
    y = 0.0
    vy = 0.0
    ay = 0.0
    return [t, x, y, z, vx, vy, vz, ax, ay, az]


def m2(t):

    r = 0.0
    T = 16.0
    w = 2.0 * np.pi / T
    if t > 4 and t < 20:
        z = 11.0 + r * (1 - np.cos(w * (t - 4.0)))
        vz = r * w * np.sin(w * (t - 4.0))
        az = r * w * w * np.cos(w * (t - 4.0))

    else:
        z = 11.0
        vz = 0.0
        az = 0.0
    y = 0.0
    vy = 0.0
    ay = 0.0
    x = 0.0
    vx = 0.0
    ax = 0.0
    return [t, x, y, z, vx, vy, vz, ax, ay, az]


n = 10000
dt = 0.01


M = np.zeros((n, 10))

for i in range(n):
    t = dt * i
    M[i, :] = m1(t)


path = "V:\\09_Mar\\01_cases\\test_case\\input\\"
np.savetxt(
    path + "dataPosicionFairlead.dat", M, fmt="%.10e", header=str(n), comments=""
)
