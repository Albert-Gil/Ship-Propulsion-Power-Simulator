"""
Holtrop–Mennen Ship Resistance and Power Prediction

Author: Albert Gil Esmendia

Description:
    This script computes total ship resistance and required propulsion
    power using the Holtrop–Mennen method for a large cruise ship
    (Royal Caribbean's Oasis-class reference geometry).

References:
    Holtrop, J., & Mennen, G.G.J. (1982). An approximate power prediction method.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ------------------------------------------------------------------------
# Physical Constants
# ------------------------------------------------------------------------
nu = 1.19e-6      # Water kinematic viscosity [m^2/s]
rho = 1025        # Water density [kg/m^3]
g = 9.81          # Gravity acceleration [m/s^2]

# ------------------------------------------------------------------------
# Ship Data (Royal Caribbean – Oasis of the Seas, reference values)
# ------------------------------------------------------------------------
v_kn = np.arange(0.005, 30.005, 0.005)  # Ship speed [knots]
v = v_kn * 0.514444                      # Ship speed [m/s]

Loa = 361.6       # Length overall [m]
Lbp = 330         # Length between perpendiculars [m]
B = 47            # Beam at waterline [m]
T = 9.1           # Draught [m]
Tf = 9.1          # Forward draught [m]
vol = 110000      # Displacement volume [m^3]

Abt = 30.6        # Bulb transverse area [m^2]
Sapp = 1000       # Wetted surface of appendages [m^2]
Arf = 2.5         # Appendage resistance factor (1+k2)
hb = 5.2          # Bow draft [m]
At = 0            # Transom area [m^2]

Cm = 0.97         # Midship coefficient
Cwp = 0.91        # Waterplane coefficient
lcb = 0           # Longitudinal center of buoyancy [% Lbp]
Cstern = 10       # Stern shape factor (0 = normal, 10 = U-shape)

propPerf = 0.73   # Propulsive efficiency
mechPerf = 0.99   # Mechanical efficiency

# ------------------------------------------------------------------------
# Derived Hull Parameters
# ------------------------------------------------------------------------
Cb = vol / (Lbp * B * T)
Cp = Cb / Cm

S = (Lbp * (2*T + B) * np.sqrt(Cm) *
     (0.453 + 0.4425*Cb - 0.2862*Cm
      - 0.003467*B/T + 0.3696*Cwp)
     + 2.38 * Abt / Cb)

Fn = v / np.sqrt(g * Lbp)
Re = v * Lbp / nu

Lr = (1 - Cp + 0.06*Cp*lcb/(4*Cp - 1)) * Lbp

ie = (1 + 89 * np.exp(
    -(Lbp/B)**0.80856
    * (1 - Cwp)**0.30484
    * (1 - Cp - 0.0225*lcb)**0.6367
    * (Lr/B)**0.34574
    * (100*vol/Lbp**3)**0.16302))

# ------------------------------------------------------------------------
# Frictional Resistance (ITTC-1957)
# ------------------------------------------------------------------------
Cf = 0.075 / (np.log10(Re) - 2)**2
Rf = 0.5 * rho * v**2 * S * Cf

# ------------------------------------------------------------------------
# Bulbous Bow Resistance
# ------------------------------------------------------------------------
Pb = 0.56 * np.sqrt(Abt) / (T - 1.5*hb)
Fni = v / np.sqrt(g*(T - hb - 0.25*np.sqrt(Abt)) + 0.15*v**2)

Rb = (0.11 * np.exp(-3 / Pb**2) * Fni**3 * Abt**(3/2) * rho * g /
      (1 + Fni**2))

# ------------------------------------------------------------------------
# Appendage Resistance
# ------------------------------------------------------------------------
Rapp = 0.5 * rho * v**2 * Sapp * Arf * Cf

# ------------------------------------------------------------------------
# Transom Resistance (not active for this hull)
# ------------------------------------------------------------------------
Rtr = np.zeros_like(v)

# ------------------------------------------------------------------------
# Wave Resistance (Holtrop–Mennen)
# ------------------------------------------------------------------------
c7 = np.zeros_like(v)
ratioBL = B / Lbp

if ratioBL <= 0.11:
    c7[:] = 0.229577 * ratioBL**(1/3)
elif ratioBL <= 0.25:
    c7[:] = ratioBL
else:
    c7[:] = 0.5 - 0.0625 / ratioBL

c1 = 2223105 * c7**3.78613 * (T/B)**1.07961 * (90 - ie)**(-1.37565)
c3 = 0.56 * Abt**1.5 / (B*T*(0.31*np.sqrt(Abt) + Tf - hb))
c2 = np.exp(-1.89 * np.sqrt(c3))
c5 = 1 - 0.8 * At / (B*T*Cm)

lambda_val = ((Lbp/B >= 12) * (1.446*Cp - 0.03*Lbp/B) +
              (Lbp/B < 12) * (1.446*Cp - 0.36))

if Cp <= 0.8:
    c16 = 8.07981*Cp - 13.8673*Cp**2 + 6.984388*Cp**3
else:
    c16 = 1.73014 - 0.7067*Cp

m1 = (0.0140407*Lbp/T - 1.75254*vol**(1/3)/Lbp
      - 4.79323*B/Lbp - c16)

c15 = -1.69385
m2 = c15 * Cp**2 * np.exp(-0.1 / Fn**2)
d = -0.9

c17 = 6919.3 * Cm**(-1.3346) * (vol/Lbp**3)**2.00977 * (Lbp/B - 2)**1.40692
m3 = -7.2035 * (B/Lbp)**0.326869 * (T/B)**0.605375
m4 = c15 * 0.4 * np.exp(-0.034 / Fn**3.29)

Rw04 = (c1 * c2 * c5 * vol * rho * g *
        np.exp(m1 * Fn**d + m4 * np.cos(lambda_val / Fn**2)))

Rw055 = (c17 * c2 * c5 * vol * rho * g *
         np.exp(m3 * Fn**d + m4 * np.cos(lambda_val / Fn**2)))

Rw = np.zeros_like(Fn)
Rw[Fn < 0.4] = Rw04[Fn < 0.4]
mask = (Fn >= 0.4) & (Fn <= 0.55)
Rw[mask] = (Rw04[mask] + (10*Fn[mask] - 4) *
            (Rw055[mask] - Rw04[mask]) / 1.5)
Rw[Fn > 0.55] = Rw055[Fn > 0.55]

# ------------------------------------------------------------------------
# Correlation Resistance
# ------------------------------------------------------------------------
c4 = min(Tf/Lbp, 0.04)
Ca = (0.006*(Lbp + 100)**(-0.16) - 0.00205
      + 0.003*np.sqrt(Lbp/7.5)*Cb**4*c2*(0.04 - c4))

Ra = 0.5 * rho * v**2 * S * Ca

# ------------------------------------------------------------------------
# Total Resistance and Power
# ------------------------------------------------------------------------
HullFactor = (0.93
              + 0.487118*(1 + 0.011*Cstern)
              * (B/Lbp)**1.06806
              * (T/Lbp)**0.46106
              * (Lbp/Lr)**0.121563
              * (Lbp**3/vol)**0.36486
              * (1 - Cp)**(-0.604247))

Rt = Rf * HullFactor + Rapp + Rw + Rb + Rtr + Ra
Rt[np.isnan(Rt)] = 0

Pe = Rt * v
Ps = Pe / (propPerf * mechPerf)

# ------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------
figPath = "./figures"
if not os.path.exists(figPath):
    os.makedirs(figPath)

# Style
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['font.size'] = 16

# Resistance plot
plt.figure(figsize=(10, 7))
plt.plot(v_kn, Rt/1e3)
plt.xlabel('Speed [kn]')
plt.ylabel('Resistance [kN]')
plt.xlim(0, 30)
plt.ylim(0, 10500)
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(figPath, "Resistance.png"))

# Power plot
plt.figure(figsize=(10, 7))
plt.plot(v_kn, Pe/1e6, label='Effective power')
plt.plot(v_kn, Ps/1e6, label='Shaft power')
plt.xlabel('Speed [kn]')
plt.ylabel('Power [MW]')
plt.xlim(0, 30)
plt.ylim(0, 220)
plt.legend(loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(figPath, "Power.png"))

plt.show()