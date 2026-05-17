"""
H-S (Mollier) Diagram using IAPWS-95 Python bindings

1. Calculating Isotherm lines isot(0.0, 800.0) °C
2. Calculating Isobar lines isop(611.657e-6, 100.0) MPa
3. Calculating saturation lines x=0, x=1
4. Calculating isoquality lines x(0.1, 0.9)

Author: Generated from seuif97 reference
"""
import numpy as np
import matplotlib.pyplot as plt
from iapws95 import pt2h, pt2s, tx2h, tx2s

xAxis = "s"
yAxis = "h"
title = {"h": "h, kJ/kg", "s": "s, kJ/(kg·K)"}

plt.title(f"{yAxis}-{xAxis} Diagram (IAPWS-95)")
plt.xlabel(title[xAxis])
plt.ylabel(title[yAxis])
plt.xlim(0, 12.5)
plt.ylim(0, 4300)
plt.grid(True, linestyle='--', alpha=0.5)

Pt = 611.657e-6

isot = np.array([0, 50, 100, 200, 300, 400, 500, 600, 700, 800])
isop = np.array([Pt, 0.001, 0.01, 0.1, 1, 10, 20, 50, 100])

for t in isot:
    h = np.array([pt2h(p, t) for p in isop])
    s = np.array([pt2s(p, t) for p in isop])
    plt.plot(s, h, 'g', lw=0.5)

for p in isop:
    h = np.array([pt2h(p, t) for t in isot])
    s = np.array([pt2s(p, t) for t in isot])
    plt.plot(s, h, 'b', lw=0.5)

tc = 647.096 - 273.15
T = np.linspace(0.1, tc, 100)

for x in np.array([0, 1.0]):
    h = np.array([tx2h(t, x) for t in T])
    s = np.array([tx2s(t, x) for t in T])
    plt.plot(s, h, 'r', lw=1.0)

isox = np.linspace(0.1, 0.9, 9)
for x in isox:
    h = np.array([tx2h(t, x) for t in T])
    s = np.array([tx2s(t, x) for t in T])
    plt.plot(s, h, 'r--', lw=0.5)

plt.tight_layout()
plt.show()
