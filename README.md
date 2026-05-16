# IAPWS-95 Thermodynamic Properties Library (Rust Implementation)

## Project Overview

This project is a Rust implementation of the [IAPWS-95](https://iapws.org/readme/iapws-r1/) (International Association for the Properties of Water and Steam Formulation 1995) standard for calculating thermodynamic properties of water and steam. The formulation was published in 1995 and revised in 2018, providing accurate calculation capabilities for thermodynamic properties of water and steam over a wide range of temperatures and pressures.

### Valid Range

| Parameter | Range |
|------|------|
| Temperature (T) | 273.16 K to 1273 K (0°C to 1000°C) |
| Pressure (p) | Up to 1000 MPa (extended range available: 100000 MPa) |
| Density (ρ) | Based on critical density ρc = 322 kg/m³ |

### Calculated Properties

- **Basic state quantities**: Temperature, density, pressure
- **Thermodynamic functions**: Internal energy (u), enthalpy (h), entropy (s)
- **Specific heat capacities**: Constant-volume specific heat (cv), constant-pressure specific heat (cp)
- **Derived properties**: Speed of sound (w), Joule-Thomson coefficient (μ), isothermal compressibility (κ)

## Quick Start

### Requirements

- Rust toolchain (edition 2021)
- Cargo build system

### Installation and Build

```bash
# Clone the project
cd rust

# Build the project
cargo build

# Run tests
cargo test

# Generate documentation
cargo doc --open
```

### Usage Example

As a library:

## Project Structure

```
rust/
├── Cargo.toml              # Project configuration
├── src/
│   ├── lib.rs              # Library entry point, exports all public modules
│   ├── iapws95.rs          # Main module: reference constants, data structures, API functions
│   ├── iapws95_ideal.rs    # Ideal gas part implementation (φ°)
│   └── iapws95_residual.rs # Residual part implementation (φʳ)
└── tests/
    ├── td_free_energy.rs   # Helmholtz free energy calculation verification test
    └── td_test.rs          # T-d-p equation of state test
```

### Module Description

#### `iapws95` - Main Module

Provides reference constants, valid range definitions, input/output data structures, and main API functions.

**Reference Constants**:
- `IAPWS95_TCRIT = 647.096 K` - Critical temperature
- `IAPWS95_RHOCRIT = 322.0 kg/m³` - Critical density
- `IAPWS95_R = 0.46151805 kJ/(kg·K)` - Specific gas constant


**Main API Functions**:

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
| `d2phi_residual_ddelta2(delta, tau)` | Second derivative ∂²φʳ/∂δ² (to be completed) |
| `dphi_residual_dtau(delta, tau)` | First derivative ∂φʳ/∂τ |
| `d2phi_residual_dtau2(delta, tau)` | Second derivative ∂²φʳ/∂τ² |
| `d2phi_residual_ddelta_dtau(delta, tau)` | Mixed derivative ∂²φʳ/∂δ∂τ |

---

## API Functions

| Function | Parameters | Description |
|------|------|------|
| `calc_pressure` | (T: f64, rho: f64) -> f64 | Calculate pressure (MPa) |
| `calc_internal_energy` | (T: f64, rho: f64) -> f64 | Calculate internal energy (kJ/kg) |
| `calc_enthalpy` | (T: f64, rho: f64, p: f64, u: f64) -> f64 | Calculate enthalpy (kJ/kg) |
| `calc_entropy` | (T: f64, phi_o: f64, phi_r: f64, dphi_dtau: f64) -> f64 | Calculate entropy (kJ/(kg·K)) |
| `calc_cv` | (T: f64, d2phi_dtau2: f64) -> f64 | Calculate constant-volume specific heat (kJ/(kg·K)) |
| `calc_cp` | (T: f64, rho: f64, dphi_ddelta: f64, d2phi_ddelta2: f64) -> f64 | Calculate constant-pressure specific heat (kJ/(kg·K)) |
| `calc_speed_of_sound` | (rho: f64, T: f64, cp: f64, cv: f64) -> f64 | Calculate speed of sound (m/s) |

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
| Internal energy u | RT·τ·(φ° + φʳ + τ·∂φ/∂τ) |
| Entropy s | R·(φ° + φʳ - τ·∂φ/∂τ) |
| Enthalpy h | u + p/ρ × 1000 |
| Constant-volume specific heat cv | R·(-τ²·∂²φ/∂τ²) |
| Constant-pressure specific heat cp | cv + R·(∂p/∂T)²ᵣₒₕₑ / (∂p/∂ρ)ₜ |
| Speed of sound w | √(cp/cv · RT × 1000) |

### Numerical Accuracy

According to the IAPWS-95 official documentation, the formulation provides in the following ranges:
- **Single-phase region**: Uncertainty < 0.01% (pressure), < 0.02% (enthalpy)
- **Saturation line**: Uncertainty < 0.05% (saturation pressure), < 0.1% (saturation density)

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
| `td_test.rs` | Verifies T-d-p equation of state calculation accuracy (based on Table 5 Region 1 reference data) | `assert_approx_eq` |

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

**T-d-p Equation of State Verification** (`td_test.rs`):
```rust
// Table 5, Page 9: Region 1 reference data (T=300K)
let Td_data = [
    propD { T: 300.0, d: 0.995660e3, p: 0.992418352e-1, cv: 0.417301218e1, w: 0.150773921e4, s: 0.393062643 },
    // ... more data points
];

for i in 0..11 {
    assert_approx_eq!(Tv_data[i].p, calc_pressure(Tv_data[i].T, Tv_data[i].v), 1.0e-6);
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

### To-Be-Completed Features

The following features need further completion. Contributions are welcome:

| Feature | Description | Priority |
|------|------|------|
| `iapws95_from_TP` | Calculate properties from temperature and pressure (requires iterative density solver) | High |
| `iapws95_saturation` | Saturation state calculation | High |
| `d2phi_residual_ddelta2` | Second derivative of residual part ∂²φʳ/∂δ² (currently a placeholder) | Medium |

#### `iapws95_from_TP` Implementation Suggestions

Requires developing an iterative density solver. Recommended approach:
1. **Initial guess**: Use ideal gas equation or Antoine equation
2. **Iteration method**: Newton-Raphson or Brent method
3. **Convergence criterion**: |p_calc - p_target| < 1e-6 MPa

#### `iapws95_saturation` Implementation Suggestions

Requires developing a phase equilibrium solver:
1. **Saturation pressure calculation**: Use Wagner equation or iterative solution
2. **Saturation density calculation**: Based on chemical potential equality condition μ'(T,p) = μ''(T,p)
3. **Range limitation**: 273.16 K to critical temperature 647.096 K

### Performance Suggestions

1. **Avoid redundant calculations**: Reduced variables δ and τ can be reused across multiple calculations
2. **Use batch processing**: For large datasets, consider calling calculation functions directly in loops

### Contribution Guidelines

1. Fork the project
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Create a Pull Request

### Code Standards

- All public APIs must have documentation comments
- New functions need test cases
- Maintain backward compatibility
- Follow IAPWS-95 official formula notation conventions

## Dependencies

### Runtime Dependencies

This project does not depend on any external crates, using only the Rust standard library:

```toml
[dependencies]
# No external dependencies
```

This ensures:
- Minimal build size
- Fastest compilation time
- Minimal security vulnerability risk
- Cross-platform compatibility

### Development Dependencies

Testing uses the following dev-dependency:

```toml
[dev-dependencies]
assert_approx_eq = "1.1.0"  # Floating-point approximate comparison macros
```

## License

This project follows the licensing requirements of the IAPWS-95 standard. Please refer to the [IAPWS official website](https://iapws.org/) for the latest licensing information.

## References

1. Wagner, W., & Pruß, A. (2002). The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use. *Journal of Physical and Chemical Reference Data*, 31(2), 387-535.

2. IAPWS (2018). Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance. [IAPWS R6-95(2018)](https://iapws.org/readme/iapws-r1/)

---

## Contact

For questions or suggestions, please contact through:
- Submit an Issue
- Create a Pull Request

