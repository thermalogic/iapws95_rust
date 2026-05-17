# IAPWS-95 Thermodynamic Properties Library (Rust Implementation)

## Project Overview

This project is a Rust implementation of the [IAPWS-95](https://iapws.org/readme/iapws-r1/) (International Association for the Properties of Water and Steam Formulation 1995) standard for calculating thermodynamic properties of water and steam. The formulation was published in 1995 and revised in 2018, providing accurate calculation capabilities for thermodynamic properties of water and steam over a wide range of temperatures and pressures.

### Valid Range

| Parameter | Range |
|------|------|
| Temperature (T) | 273.16 K to 1273 K (0°C to 1000°C) |
| Pressure (p) | Up to 1000 MPa (extended range: 100000 MPa) |
| Density (ρ) | Based on critical density ρc = 322 kg/m³ |

## Project Status

> ⚠️ **This project is under active development.** The following sections describe the current implementation status.

### ✅ Implemented Features

| Category | Feature | Status |
|------|------|------|
| **Core Formulation** | Ideal gas part φ°(δ,τ) | ✅ Complete |
| | Residual part φʳ(δ,τ) - Polynomial terms | ✅ Complete |
| | Residual part φʳ(δ,τ) - Exponential terms | ✅ Complete |
| | Residual part φʳ(δ,τ) - Gaussian terms | ✅ Complete |
| | Residual part φʳ(δ,τ) - Non-analytic terms | ✅ Complete |
| **Derivatives** | ∂φ/∂δ, ∂²φ/∂δ² | ✅ Complete |
| | ∂φ/∂τ, ∂²φ/∂τ² | ✅ Complete |
| | ∂²φ/∂δ∂τ | ✅ Complete |
| **Properties** | Pressure (p) | ✅ Verified |
| | Constant-volume specific heat (cv) | ✅ Verified |
| | Speed of sound (w) | ✅ Verified (error < 0.2 m/s near critical point) |
| | Entropy (s) | ✅ Verified |
| | Internal energy (u) | ✅ Implemented |
| | Enthalpy (h) | ✅ Implemented |
| | Constant-pressure specific heat (cp) | ✅ Implemented |
| **Saturation** | Saturation properties module | ✅ Complete |
| | Saturation pressure pₛ(T) | ✅ Complete |
| | Saturation densities ρ'(T), ρ''(T) | ✅ Complete |
| | Saturation enthalpy h'(T), h''(T) | ✅ Complete |
| | Saturation entropy s'(T), s''(T) | ✅ Complete |
| **Testing** | Helmholtz free energy verification | ✅ Complete |
| | T-ρ-p property verification (11 test points) | ✅ Complete |

### 🚧 Planned Features

| Category | Feature | Priority |
|------|------|------|
| **Properties** | Isentropic exponent | 🔜 |
| | Joule-Thomson coefficient | 🔜 |
| | Thermal expansion coefficient | 🔜 |
| | Isothermal compressibility | 🔜 |
| **Validation** | Extended test coverage (all IAPWS-95 tables) | 🔜 |
| | Benchmark against reference implementations | 🔜 |

Legend: 🔜 = Near term, 🔮 = Future consideration

## Quick Start

### Requirements

**For Rust usage**:
- Rust toolchain (edition 2021)
- Cargo build system

**For Python bindings**:
- Python 3.8+
- maturin >= 1.0
- pip package manager

**For C/C++ bindings**:
- C compiler (gcc, clang, or MSVC)
- Rust toolchain for building the library

### Installation and Build

#### Building the Rust Library

```bash
# Clone the project
cd iapws95_rust

# Build the project (debug mode)
cargo build

# Build release version (optimized)
cargo build --release

# Run tests
cargo test

# Generate documentation
cargo doc --open
```

#### Building Python Bindings

```bash
# Create and activate virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Build and install the package
maturin develop --features python

# Or build wheel package for distribution
maturin build --features python
pip install target/wheels/iapws95-*.whl
```

#### Building C/C++ Bindings

```bash
# Build release version (produces shared library)
cargo build --release

# The library will be available at:
# Linux:   target/release/libiapws95.so
# Windows: target/release/iapws95.dll + iapws95.lib
# macOS:   target/release/libiapws95.dylib
```

### Usage Example

#### Rust (Direct Computation)

```rust
use iapws95::iapws95::*;

// Calculate properties at T=500K, rho=838.025 kg/m³
let T = 500.0;
let rho = 838.025;

let p = calc_pressure(T, rho);
let cv = calc_cv(T, rho);
let w = calc_speed_of_sound(T, rho);
let s = calc_entropy(T, rho);
```

#### Python (Direct Computation)

```python
from iapws95 import tr2h, tr2s, tr2cp

# Calculate properties at T=500°C, ρ=838.025 kg/m³
t = 500.0
rho = 838.025

h = tr2h(t, rho)      # Enthalpy: kJ/kg
s = tr2s(t, rho)      # Entropy: kJ/(kg·K)
cp = tr2cp(t, rho)    # Cp: kJ/(kg·K)
```

#### C/C++ (Direct Computation)

```c
#include "iapws95.h"

// Calculate properties at T=500°C, ρ=838.025 kg/m³
double t = 500.0;
double rho = 838.025;

double h = iapws95_tr2h(t, rho);    // Enthalpy: kJ/kg
double s = iapws95_tr2s(t, rho);    // Entropy: kJ/(kg·K)
double cp = iapws95_tr2cp(t, rho);  // Cp: kJ/(kg·K)
```

## Project Structure

```
iapws95_rust/
├── Cargo.toml              # Project configuration (includes pyo3 for Python bindings)
├── pyproject.toml          # maturin build configuration for Python package
├── src/
│   ├── lib.rs                    # Library entry point, exports all public modules
│   ├── iapws95.h                 # C FFI header file
│   ├── iapws95.rs                # Main module: reference constants, data structures, API functions
│   ├── iapws95_ideal.rs          # Ideal gas part implementation (φ°)
│   ├── iapws95_residual.rs       # Residual part implementation (φʳ)
│   ├── iapws95_saturation.rs     # Saturation properties calculation module
│   ├── py_iapws95.rs             # Python bindings via PyO3
│   └── c_iapws95.rs              # C FFI bindings
├── examples/
│   └── basic_usage.rs            # Example: single-phase and saturation properties
├── demo/
│   ├── iapws95_usage.py          # Python example: Basic property calculation
│   ├── c_example.c               # C example: Using C FFI bindings
│   └── Makefile                  # Build script for C examples
└── tests/
    ├── td_free_energy.rs         # Helmholtz free energy calculation verification test
    ├── td_test.rs                # T-ρ-p equation of state test (Table 7)
    └── T_saturation_table8.rs    # Saturation properties verification (Table 8)
```

### Module Description

#### `iapws95` - Main Module

Provides reference constants, valid range definitions, and main API functions.

**Reference Constants**:
- `IAPWS95_TCRIT = 647.096 K` - Critical temperature
- `IAPWS95_RHOCRIT = 322.0 kg/m³` - Critical density
- `IAPWS95_PCRIT = 22.064 MPa` - Critical pressure
- `IAPWS95_R = 0.46151805 kJ/(kg·K)` - Specific gas constant

**Helper Functions**:

| Function | Description |
|------|------|
| `reduced_density(rho)` | Calculate reduced density δ = ρ/ρc |
| `inv_reduced_temp(T)` | Calculate inverse reduced temperature τ = Tc/T |
| `iapws95_in_range(T, p)` | Check if state is within valid range |

**Main API Functions**:

| Function | Description |
|------|------|
| `calc_pressure(T, rho)` | Calculate pressure (MPa) |
| `calc_internal_energy(T, rho)` | Calculate internal energy (kJ/kg) |
| `calc_enthalpy(T, rho)` | Calculate enthalpy (kJ/kg) |
| `calc_entropy(T, rho)` | Calculate entropy (kJ/(kg·K)) |
| `calc_cv(T, rho)` | Calculate constant-volume specific heat (kJ/(kg·K)) |
| `calc_cp(T, rho)` | Calculate constant-pressure specific heat (kJ/(kg·K)) |
| `calc_speed_of_sound(T, rho)` | Calculate speed of sound (m/s) |

#### `iapws95_ideal` - Ideal Gas Part

Implements the ideal gas part (φ°) of the dimensionless Helmholtz free energy, based on IAPWS-95 Equation 5 and Tables 1 and 4.

**Calculation Formula**:

```
φ°(δ,τ) = ln(δ) + n₁ + n₂τ + n₃ln(τ) + Σᵢ₌₄⁸ nᵢln[1-exp(-γᵢτ)]
```

Where:
- `δ = ρ/ρc` - Reduced density
- `τ = Tc/T` - Inverse reduced temperature

**Available Functions**:

| Function | Description | Formula |
|------|------|------|
| `phi_ideal(delta, tau)` | Helmholtz free energy | φ°(δ,τ) |
| `dphi_ideal_ddelta(delta)` | First derivative ∂φ°/∂δ | 1/δ |
| `d2phi_ideal_ddelta2(delta)` | Second derivative ∂²φ°/∂δ² | -1/δ² |
| `dphi_ideal_dtau(tau)` | First derivative ∂φ°/∂τ | n₂ + n₃/τ + Σ... |
| `d2phi_ideal_dtau2(tau)` | Second derivative ∂²φ°/∂τ² | -n₃/τ² - Σ... |
| `d2phi_ideal_dtaudelta(delta, tau)` | Mixed derivative ∂²φ°/∂τ∂δ | 0 |

#### `iapws95_residual` - Residual Part

Implements the residual part (φʳ) of the dimensionless Helmholtz free energy, based on coefficients from IAPWS-95 Table 5.

**Calculation Formula**:

```
φʳ(δ,τ) = Σᵢ nᵢδᵈⁱτᵗⁱ                                    [Polynomial terms, i=1-7]
        + Σᵢ nᵢδᵈⁱτᵗⁱexp(-δᶜⁱ)                          [Exponential terms, i=8-51]
        + Σᵢ nᵢδᵈⁱτᵗⁱexp[-αᵢ(δ-εᵢ)²-βᵢ(τ-γᵢ)²]          [Gaussian terms, i=52-54]
        + Σᵢ nᵢΔᵇⁱδF(δ,τ)                               [Non-analytic terms, i=55-56]
```

**Term Classification**:

| Term Type | Index Range | c value | Count |
|------|------|------|------|
| Polynomial | i=1-7 | - | 7 |
| Exponential (c=1) | i=8-22 | 1 | 15 |
| Exponential (c=2) | i=23-42 | 2 | 20 |
| Exponential (c=3,4,6) | i=43-51 | 3,4,6 | 9 |
| Gaussian | i=52-54 | - | 3 |
| Non-analytic | i=55-56 | - | 2 |

**Available Functions**:

| Function | Description |
|------|------|
| `phi_residual(delta, tau)` | Helmholtz free energy φʳ |
| `dphi_residual_ddelta(delta, tau)` | First derivative ∂φʳ/∂δ |
| `d2phi_residual_ddelta2(delta, tau)` | Second derivative ∂²φʳ/∂δ² |
| `dphi_residual_dtau(delta, tau)` | First derivative ∂φʳ/∂τ |
| `d2phi_residual_dtau2(delta, tau)` | Second derivative ∂²φʳ/∂τ² |
| `d2phi_residual_ddelta_dtau(delta, tau)` | Mixed derivative ∂²φʳ/∂δ∂τ |

#### `iapws95_saturation` - Saturation Properties Module

Implements saturation properties calculation along the vapor-liquid equilibrium line based on IAPWS-95 Table 8.

**Algorithm**:

The module uses a hybrid approach combining IAPWS SR1-86 (1992) explicit equations with Newton's method for IAPWS-95 phase equilibrium:

1. **Initial guesses**: IAPWS SR1-86 (1992) explicit equations for saturated liquid and vapor densities
   - Wagner-type correlation for liquid density: ρ'/ρc = 1 + b₁τ^(1/3) + b₂τ^(2/3) + ...
   - Exponential correlation for vapor density: ln(ρ''/ρc) = c₁τ^(2/6) + c₂τ^(4/6) + ...

2. **Phase equilibrium refinement**: Newton's method solving the IAPWS-95 phase equilibrium conditions:
   ```
   F1(δL, δV) = K(δV, τ) - K(δL, τ) = 0    (equal chemical potential)
   F2(δL, δV) = J(δV, τ) - J(δL, τ) = 0    (equal pressure)
   ```
   Where:
   - `J = δ·(1 + δ·∂φʳ/∂δ)` — dimensionless pressure term
   - `K = δ·∂φʳ/∂δ + φʳ + ln(δ)` — dimensionless chemical potential term

3. **Property calculation**: Once δ' and δ'' are determined, all saturation properties are computed using IAPWS-95 formulas

**Phase Equilibrium Condition**:

The saturation properties satisfy the Maxwell criterion:

```
p(δ', τ) = p(δ'', τ) = p_σ     (equal pressure)
μ(δ', τ) = μ(δ'', τ)           (equal chemical potential)
```

Where:
- `δ'` - Reduced density of saturated liquid
- `δ''` - Reduced density of saturated vapor
- `p_σ` - Saturation vapor pressure
- `μ` - Chemical potential

**Available Functions**:

| Function | Description | Returns |
|------|------|------|
| `calc_saturation_properties(T)` | Calculate all saturation properties at temperature T | `Option<SaturationProperties>` |

**SaturationProperties Structure**:

| Field | Description | Unit |
|------|------|------|
| `p_sat` | Saturation vapor pressure | MPa |
| `rho_l` | Saturated liquid density | kg/m³ |
| `rho_v` | Saturated vapor density | kg/m³ |
| `h_l` | Saturated liquid specific enthalpy | kJ/kg |
| `h_v` | Saturated vapor specific enthalpy | kJ/kg |
| `s_l` | Saturated liquid specific entropy | kJ/(kg·K) |
| `s_v` | Saturated vapor specific entropy | kJ/(kg·K) |

**Valid Range**:

| Parameter | Range |
|------|------|
| Temperature | 273.16 K to 647.096 K (triple point to critical point) |

---

#### `py_iapws95` - Python Bindings Module

Provides Python bindings via PyO3, enabling thermodynamic property calculations from Python code. This module is compiled as a native extension using maturin.

**Building and Installation**:

```bash
# Build and install to current virtual environment
maturin develop --features python

# Or build wheel package
maturin build --features python
pip install target/wheels/iapws95-*.whl
```

**Function Categories**:

The Python bindings provide three categories of functions based on input parameters:

**(T,ρ) → Property Calculations (Direct Computation)**:

| Function | Description | Parameters | Returns |
|------|------|------|------|
| `tr2p(t_c, rho)` | Calculate pressure at given temperature and density | t: °C, ρ: kg/m³ | p: MPa |
| `tr2u(t_c, rho)` | Calculate internal energy at given temperature and density | t: °C, ρ: kg/m³ | u: kJ/kg |
| `tr2h(t_c, rho)` | Calculate enthalpy at given temperature and density | t: °C, ρ: kg/m³ | h: kJ/kg |
| `tr2s(t_c, rho)` | Calculate entropy at given temperature and density | t: °C, ρ: kg/m³ | s: kJ/(kg·K) |
| `tr2cv(t_c, rho)` | Calculate constant-volume specific heat at given T and ρ | t: °C, ρ: kg/m³ | cv: kJ/(kg·K) |
| `tr2cp(t_c, rho)` | Calculate constant-pressure specific heat at given T and ρ | t: °C, ρ: kg/m³ | cp: kJ/(kg·K) |
| `tr2w(t_c, rho)` | Calculate speed of sound at given temperature and density | t: °C, ρ: kg/m³ | w: m/s |

**Saturation Properties**:

| Function | Description | Parameters | Returns |
|------|------|------|------|
| `saturation_properties(t)` | Calculate all saturation properties at temperature | t: °C | (p_sat, ρ', ρ'', h', h'', s', s'') |

**Usage Example**:

```python
from iapws95 import tr2h, tr2s, tr2cp, tr2cv, tr2w

# Direct calculation at T=500°C, ρ=838.025 kg/m³ (recommended approach)
t = 500.0
rho = 838.025

h = tr2h(t, rho)      # Enthalpy: kJ/kg
s = tr2s(t, rho)      # Entropy: kJ/(kg·K)
cp = tr2cp(t, rho)    # Cp: kJ/(kg·K)
cv = tr2cv(t, rho)    # Cv: kJ/(kg·K)
w = tr2w(t, rho)      # Speed of sound: m/s

# Saturation properties at T=100°C
p_sat, rho_l, rho_v, h_l, h_v, s_l, s_v = saturation_properties(100.0)

# Two-phase mixture at T=200°C, x=0.5 (quality)
h_mix = tx2h(200.0, 0.5)
s_mix = tx2s(200.0, 0.5)
```

**Dependencies**:

Python bindings require:
- Python 3.8+
- matplotlib (for plotting examples)
- numpy (for plotting examples)

**Known Limitations**:

The `(p,T)` → property functions (`pt2h()`, `pt2s()`) use numerical inversion to solve for density from pressure and temperature. This may fail in certain conditions:
- Very low pressures (< 0.001 MPa) near triple point temperature
- Near saturation boundary where multiple density solutions exist
- Extreme superheated vapor conditions

**Recommendation**: For reliable calculations, use `(T,ρ)` → property functions (`tr2p()`, `tr2h()`, etc.) which provide direct computation without numerical inversion.

---

#### `c_iapws95` - C Language Bindings Module

Provides a C-compatible interface via FFI, enabling thermodynamic property calculations from C/C++ code. The library is compiled as a shared/static library.

**Calling Convention**:

All exported functions use the standard **cdecl (C declaration)** calling convention:
- Parameters are passed on the stack from right to left
- The caller is responsible for cleaning up the stack after the call
- Compatible with all major C/C++ compilers and platforms
- No special compiler directives or pragmas required when linking

This ensures maximum portability across:
- **Windows**: MSVC, MinGW-w64, Cygwin
- **Linux**: gcc, clang
- **macOS**: clang (Xcode)

**Building the Library**:

```bash
# Build release version (produces .so/.dll/.dylib)
cargo build --release
```

This produces:
- **Linux**: `target/release/libiapws95.so`
- **Windows**: `target/release/iapws95.dll` + `target/release/iapws95.lib`
- **macOS**: `target/release/libiapws95.dylib`

**Function Categories**:

The C bindings provide three categories of functions based on input parameters:

**(T,ρ) → Property Calculations (Direct Computation)**:

| Function | Description | Parameters | Returns |
|------|------|------|------|
| `iapws95_tr2p(t_c, rho)` | Calculate pressure at given temperature and density | t: °C, ρ: kg/m³ | p: MPa |
| `iapws95_tr2u(t_c, rho)` | Calculate internal energy at given T and ρ | t: °C, ρ: kg/m³ | u: kJ/kg |
| `iapws95_tr2h(t_c, rho)` | Calculate enthalpy at given temperature and density | t: °C, ρ: kg/m³ | h: kJ/kg |
| `iapws95_tr2s(t_c, rho)` | Calculate entropy at given T and ρ | t: °C, ρ: kg/m³ | s: kJ/(kg·K) |
| `iapws95_tr2cv(t_c, rho)` | Calculate Cv at given temperature and density | t: °C, ρ: kg/m³ | cv: kJ/(kg·K) |
| `iapws95_tr2cp(t_c, rho)` | Calculate Cp at given T and ρ | t: °C, ρ: kg/m³ | cp: kJ/(kg·K) |
| `iapws95_tr2w(t_c, rho)` | Calculate speed of sound at given T and ρ | t: °C, ρ: kg/m³ | w: m/s |

**Saturation Properties**:

| Function | Description | Parameters | Returns |
|------|------|------|------|
| `iapws95_saturation_properties(t_c, props)` | Calculate all saturation properties at temperature | t: °C, props: struct* | 0 on success, -1 on error |
| `iapws95_version()` | Get library version string | - | const char* |

**Usage Example**:

```c
#include <stdio.h>
#include "iapws95.h"

int main() {
    // Direct calculation at T=500°C, ρ=838.025 kg/m³ (recommended approach)
    double t = 500.0;
    double rho = 838.025;

    double p = iapws95_tr2p(t, rho);      // Pressure: MPa
    double u = iapws95_tr2u(t, rho);      // Internal energy: kJ/kg
    double h = iapws95_tr2h(t, rho);      // Enthalpy: kJ/kg
    double s = iapws95_tr2s(t, rho);      // Entropy: kJ/(kg·K)
    double cv = iapws95_tr2cv(t, rho);    // Cv: kJ/(kg·K)
    double cp = iapws95_tr2cp(t, rho);    // Cp: kJ/(kg·K)
    double w = iapws95_tr2w(t, rho);      // Speed of sound: m/s

    printf("Pressure: %.6f MPa\n", p);
    printf("Enthalpy: %.4f kJ/kg\n", h);
    printf("Entropy: %.6f kJ/(kg·K)\n", s);

    // Saturation properties at T=100°C
    iapws95_saturation_props_t sat;
    if (iapws95_saturation_properties(100.0, &sat) == 0) {
        printf("Saturation pressure: %.6f MPa\n", sat.p_sat);
    }
    return 0;
}
```

**Compilation**:

Linux/Mac:
```bash
gcc -I../src c_example.c -L../target/release -liapws95 -o iapws95_demo
```

Windows (MSVC):
```cmd
cl /I..\src c_example.c ..\target\release\iapws95.lib
```

**Dependencies**:
- C compiler (gcc, clang, or MSVC)
- No additional runtime dependencies

---

## Algorithm Description

### Helmholtz Free Energy Formula

The core of IAPWS-95 is the dimensionless Helmholtz free energy φ(δ,τ), which is decomposed into two parts:

```
φ(δ,τ) = φ°(δ,τ) + φʳ(δ,τ)
```

Where:
- **φ°(δ,τ)**: Ideal gas part, depends only on temperature
- **φʳ(δ,τ)**: Residual part, accounts for intermolecular interactions

#### Reduced Variables

```
δ = ρ/ρc    (reduced density)
τ = Tc/T    (inverse reduced temperature)
```

### Property Calculation

All thermodynamic properties are calculated from the Helmholtz free energy and its derivatives (based on IAPWS-95 Table 3):

| Property | Formula |
|------|------|
| Pressure p | RT·δ·(1 + δ·∂φʳ/∂δ) |
| Internal energy u | RT·τ·(∂φ/∂τ) |
| Entropy s | R·(τ·∂φ/∂τ - φ° - φʳ) |
| Enthalpy h | RT·[τ·(∂φ/∂τ) + 1 + δ·(∂φʳ/∂δ)] |
| Constant-volume specific heat cv | R·(-τ²·∂²φ/∂τ²) |
| Constant-pressure specific heat cp | cv + R·(1 + δ·∂φʳ/∂δ - δ·τ·∂²φʳ/∂δ∂τ)² / (1 + 2δ·∂φʳ/∂δ + δ²·∂²φʳ/∂δ²) |
| Speed of sound w | √(RT·[1 + 2δ·∂φʳ/∂δ + δ²·∂²φʳ/∂δ² - (1 + δ·∂φʳ/∂δ - δ·τ·∂²φʳ/∂δ∂τ)² / (τ²·∂²φ/∂τ²)]) |

### Numerical Accuracy

According to the IAPWS-95 official documentation, the formulation provides the following accuracy:
- **Single-phase region**: Uncertainty < 0.01% (pressure), < 0.02% (enthalpy)
- **Saturation line**: Uncertainty < 0.05% (saturation pressure), < 0.1% (saturation density)

## Examples

### Running Rust Examples

```bash
# Run the basic usage example
cargo run --example basic_usage

# Generate Mollier (H-S) diagram using plotters
cargo run --example mollier_diagram
```

The `mollier_diagram` example generates a high-quality H-S diagram with:
- **Isotherm lines**: 0°C to 800°C (green curves)
- **Isobar lines**: 611.657 μPa to 100 MPa (blue curves)
- **Saturation dome**: x=0 (saturated liquid) and x=1 (saturated vapor) in red
- **Isoquality lines**: x=0.1 to 0.9 (red dashed curves)

Output: `mollier_diagram.png` (1200×900 pixels)

### Running Python Examples

```bash
# Build and install the Python package
maturin develop --features python

# Run basic usage example
python demo/iapws95_usage.py

```

The Python examples demonstrate:
- **Basic property calculation**: Calculate enthalpy and entropy at given (p, T) conditions using `iapws95_usage.py`

### Running C Examples

```bash
# Build the Rust library first
cargo build --release

# Compile and run the example (Linux/Mac)
cd demo
gcc -I../src c_example.c -L../target/release -liapws95 -o iapws95_demo
./iapws95_demo

# Or use Makefile
make build
make run
```

The C example demonstrates:
- **Saturation properties**: Get all saturation properties at a given temperature
- **Multi-temperature table**: Generate a saturation properties table from triple point to critical point
- **Direct (T,ρ) calculations**: Compute pressure, internal energy, enthalpy, entropy, specific heats, and speed of sound directly without numerical inversion

**Note**: The C example includes comprehensive demonstrations of all function categories, including the recommended `(T,ρ)` → property direct computation functions.

---

## Testing

### Running Tests

```bash
# Run all tests
cargo test

# Run specific test files
cargo test --test td_free_energy
cargo test --test td_test

# Show test output
cargo test -- --nocapture
```

### Test Coverage

The current test suite includes the following test cases:

| Test File | Description | Dependencies |
|------|------|------|
| `td_free_energy.rs` | Verifies calculation accuracy of ideal gas and residual parts of Helmholtz free energy (based on Table 6 reference data) | `assert_approx_eq` |
| `td_test.rs` | Verifies T-d-p equation of state calculation accuracy (based on Table 7 reference data) | `assert_approx_eq` |
| `T_saturation_table8.rs` | Verifies saturation properties calculation accuracy (based on Table 8 reference data) | `assert_approx_eq` |

### Test Case Examples

**Helmholtz Free Energy Verification** (`td_free_energy.rs`):
```rust
// Table 6: T=500K, ρ=838.025 kg/m³
let t_test = 500.0;
let rho_test = 838.025;

let delta = reduced_density(rho_test);
let tau = inv_reduced_temp(t_test);

let phi_o_ref = 2.04797733;   // φ° reference value
let phi_r_ref = -3.42693206;  // φʳ reference value

assert_approx_eq!(phi_o_ref, phi_ideal(delta, tau), 1.0e-6);
assert_approx_eq!(phi_r_ref, phi_residual(delta, tau), 1.0e-6);
```

**T-ρ-p Equation of State Verification** (`td_test.rs`):
```rust
// Table 7: T-d-p reference data
let Td_data = [
    propD { T: 300.0, d: 0.9965560e3, p: 0.992418352e-1, cv: 0.413018112e1, w: 0.150151914e4, s: 0.393062643 },
    // ... more data points
];

for i in 0..11 {
    assert_approx_eq!(Td_data[i].p, calc_pressure(Td_data[i].T, Td_data[i].d), 1.0e-6);
    assert_approx_eq!(Td_data[i].cv, calc_cv(Td_data[i].T, Td_data[i].d), 1.0e-4);
    assert_approx_eq!(Td_data[i].w, calc_speed_of_sound(Td_data[i].T, Td_data[i].d), 1.0e-5);
    assert_approx_eq!(Td_data[i].s, calc_entropy(Td_data[i].T, Td_data[i].d), 1.0e-6);
}
```

---

## Development Guide

### Code Style

This project follows Rust community coding conventions:
- Use `rustfmt` to format code
- Add documentation comments (`///`) for public APIs
- Use `#[inline]` to optimize frequently called functions
- Leverage Rust's type system for safety

### Adding New Property Calculations

To add new thermodynamic properties, you need to:

1. **Add calculation formula in `iapws95.rs`**:
```rust
pub fn calc_new_property(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    // ... calculation based on Helmholtz free energy derivatives
}
```

### Code Standards

- All public APIs must have documentation comments
- New functions need test cases
- Maintain backward compatibility
- Follow IAPWS-95 official formula notation conventions

## Dependencies

### Runtime Dependencies (Rust Core Library)

The core Rust library has **no external dependencies**, using only the Rust standard library:

```toml
[dependencies]
# No external dependencies
```

This ensures:
- Minimal build size and fastest compilation time
- Minimal security vulnerability surface
- Maximum cross-platform compatibility
- Zero dependency tree to manage

### Python Bindings Dependencies

For building Python bindings (optional feature):

**Build-time**:
```toml
[dependencies]
pyo3 = { version = "0.24", features = ["extension-module"], optional = true }

[features]
python = ["pyo3"]
cffi = []  # C FFI bindings
all-bindings = ["python", "cffi"]  # Enable all bindings
```

Build with: `cargo build --features python` or `maturin develop --features python`

**Runtime**:
- Python 3.8+
- matplotlib >= 3.0 (for plotting examples)
- numpy >= 1.20 (for plotting examples)

Install Python dependencies via pip:
```bash
pip install matplotlib numpy
```

### C/C++ Bindings Dependencies

For building C/C++ bindings:
- **Build-time**: Rust toolchain with `cargo build --release`
- **Runtime**: No additional dependencies required
- **Compiler**: gcc, clang, or MSVC (for compiling example code)

The shared library (`libiapws95.so`, `iapws95.dll`, or `libiapws95.dylib`) is self-contained and requires no external libraries.

### Build Tools

Required build tools:
- **maturin >= 1.0**: For building Python packages from Rust code
- **Cargo**: Rust's package manager and build system (included with Rust toolchain)

Optional tools:
- **pip**: Python package installer (for installing dependencies)
- **gcc/clang/MSVC**: C/C++ compilers (for compiling example code)

### Development Dependencies

Testing and development use the following dev-dependencies:

```toml
[dev-dependencies]
assert_approx_eq = "1.1.0"  # Floating-point approximate comparison macros
plotters = { version = "0.3", features = ["all_series"] }  # For H-S diagram generation in examples
```

## Features

The library provides optional features for different binding types:

| Feature | Description | Dependencies | Use Case |
|------|------|------|------|
| `python` | Enable Python bindings via PyO3 | pyo3 >= 0.24 | Building Python extension modules |
| `cffi` | Enable C FFI bindings | None (always available) | Building C/C++ shared libraries |
| `all-bindings` | Enable all optional bindings | pyo3 + cffi | Full multi-language support |

**Usage Examples**:

```bash
# Build with Python bindings only
cargo build --features python

# Build with C FFI bindings only (default for non-test builds)
cargo build --features cffi

# Build with all bindings
cargo build --features all-bindings

# Build Python package with maturin
maturin develop --features python
```

**Default Behavior**:
- Core library: Always available, no features required
- C FFI bindings: Enabled by default in non-test builds (for backward compatibility)
- Python bindings: Must be explicitly enabled with `--features python`

## License

This project follows the licensing requirements of the IAPWS-95 standard. Please refer to the [IAPWS official website](https://iapws.org/) for the latest licensing information.

## References

1. Wagner, W., & Pruß, A. (2002). The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use. *Journal of Physical and Chemical Reference Data*, 31(2), 387-535.

2. IAPWS (2018). Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance. [IAPWS R6-95(2018)](https://iapws.org/readme/iapws-r1/)
