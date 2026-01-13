# Ship Resistance and Power Simulator

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

This repository provides **clean, vectorized MATLAB and Python implementations** of the Holtrop–Mennen empirical method for estimating ship resistance and propulsion power.  
The reference case corresponds to a **Royal Caribbean Oasis-class cruise ship**, but the scripts can be adapted to other displacement vessels by modifying the input parameters.

---

## 📌 Features

- ITTC-1957 frictional resistance formulation  
- Full Holtrop–Mennen wave resistance model  
- Bulbous bow, appendage, and correlation resistance terms  
- Effective power and shaft power estimation  
- Both MATLAB and Python implementations with identical results

---

## ⚙️ Requirements

### MATLAB
- MATLAB R2020a or newer  
- No additional toolboxes required

### Python
- Python 3.7 or newer
- NumPy
- Matplotlib

Install Python dependencies:
```bash
pip install numpy matplotlib
```

---

## ▶️ Usage

### MATLAB

1. Clone or download the repository  
2. Open MATLAB in the project directory  
3. Run:

```matlab
PropulsionResistancePower
```

### Python

1. Clone or download the repository  
2. Navigate to the project directory  
3. Run:

```bash
python PropulsionResistancePower.py
```

Both scripts will:
- Compute resistance and power over a speed range (0-30 knots)
- Generate resistance and power curves
- Save figures to the `figures/` folder

---

## 🛠️ Modifying the Ship Geometry

Edit the **Ship Data** section in either `PropulsionResistancePower.m` or `PropulsionResistancePower.py`:

**MATLAB:**
```matlab
Lbp = 330;    % Length between perpendiculars [m]
B   = 47;     % Beam [m]
T   = 9.1;    % Draught [m]
vol = 110000; % Displacement volume [m^3]
```

**Python:**
```python
Lbp = 330     # Length between perpendiculars [m]
B = 47        # Beam [m]
T = 9.1       # Draught [m]
vol = 110000  # Displacement volume [m^3]
```

All coefficients are derived automatically from these parameters.

---

## 📊 Output

Both implementations generate two plots:

1. **Resistance.png**: Total ship resistance vs. speed
2. **Power.png**: Effective power and shaft power vs. speed

Plots are saved in the `figures/` directory with consistent formatting between MATLAB and Python versions.

---

## 📖 References

- Holtrop, J., & Mennen, G. G. J. (1978). *A Statistical Power Prediction Method.*  
- Holtrop, J., & Mennen, G. G. J. (1982). *An Approximate Power Prediction Method.*   
- Kristensen, H. O., & Lützen, M. (2012). *Prediction of Resistance and Propulsion Power of Ships.*  

---

## 👤 Author

Albert Gil Esmendia
