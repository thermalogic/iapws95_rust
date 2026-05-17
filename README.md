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

### Usage Example

#### Rust

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

## Project Structure

```
iapws95_rust/
├── Cargo.toml              # Project configuration
├── src/
│   ├── lib.rs                    # Library entry point, exports all public modules
│   ├── iapws95.rs                # Main module: reference constants, data structures, API functions
│   ├── iapws95_ideal.rs          # Ideal gas part implementation (φ°)
│   ├── iapws95_residual.rs       # Residual part implementation (φʳ)
│   ├── iapws95_saturation.rs     # Saturation properties calculation module
├── examples/
│   └── basic_usage.rs            # Example: single-phase and saturation properties
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


**(t,ρ) → Property Calculations (Direct Computation)**:

| Function | Description | Parameters | Returns |
|------|------|------|------|
| `tr2p(t_c, rho)` | Calculate pressure at given temperature and density | t: °C, ρ: kg/m³ | p: MPa |
| `tr2u(t_c, rho)` | Calculate internal energy at given temperature and density | t: °C, ρ: kg/m³ | u: kJ/kg |
| `tr2h(t_c, rho)` | Calculate enthalpy at given temperature and density | t: °C, ρ: kg/m³ | h: kJ/kg |
| `tr2s(t_c, rho)` | Calculate entropy at given temperature and density | t: °C, ρ: kg/m³ | s: kJ/(kg·K) |
| `tr2cv(t_c, rho)` | Calculate constant-volume specific heat at given T and ρ | t: °C, ρ: kg/m³ | cv: kJ/(kg·K) |
| `tr2cp(t_c, rho)` | Calculate constant-pressure specific heat at given T and ρ | t: °C, ρ: kg/m³ | cp: kJ/(kg·K) |
| `tr2w(t_c, rho)` | Calculate speed of sound at given temperature and density | t: °C, ρ: kg/m³ | w: m/s |

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
| `sat_t(t）` | Calculate all saturation properties at temperature t,°C | `Option<SaturationProperties>` |

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
```

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

### Development Dependencies

Testing and development use the following dev-dependencies:

```toml
[dev-dependencies]
assert_approx_eq = "1.1.0"  # Floating-point approximate comparison macros
```

## License

This project follows the licensing requirements of the IAPWS-95 standard. Please refer to the [IAPWS official website](https://iapws.org/) for the latest licensing information.

## References

1. Wagner, W., & Pruß, A. (2002). The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use. *Journal of Physical and Chemical Reference Data*, 31(2), 387-535.

2. IAPWS (2018). Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance. [IAPWS R6-95(2018)](https://iapws.org/readme/iapws-r1/)
