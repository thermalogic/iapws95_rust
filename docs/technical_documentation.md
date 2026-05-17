# IAPWS-95 Thermodynamic Properties Library — Technical Documentation

**Version:** 0.1.0  
**Language:** Rust (Edition 2021)  
**Standard:** IAPWS-95 (Revised 2018)

---

## Table of Contents

1. [Overview](#1-overview)
2. [Theoretical Foundation](#2-theoretical-foundation)
3. [Project Structure](#3-project-structure)
4. [Module Reference](#4-module-reference)
5. [API Reference](#5-api-reference)
6. [Saturation Properties Algorithm](#6-saturation-properties-algorithm)
7. [Test Suite](#7-test-suite)
8. [Build and Usage](#8-build-and-usage)

---

## 1. Overview

This library implements the **IAPWS-95** formulation — the International Association for the Properties of Water and Steam standard for thermodynamic properties of ordinary water substance. It provides accurate calculations of pressure, enthalpy, entropy, specific heats, speed of sound, and saturation properties for water and steam.

### Key Features

- Zero external runtime dependencies (Rust standard library only)
- Full IAPWS-95 Helmholtz free energy formulation (56 terms across 4 categories)
- Single-phase property calculations in the (T, ρ) domain
- Saturation properties with phase equilibrium solver
- Verified against IAPWS reference tables (Tables 6, 7, 8)

### Valid Range

| Parameter | Range |
|-----------|-------|
| Temperature T | 273.16 K – 1273 K (0°C – 1000°C) |
| Pressure p | Up to 1000 MPa (extended: 100,000 MPa) |

### Reference Constants

| Constant | Symbol | Value | Unit |
|----------|--------|-------|------|
| Critical temperature | Tc | 647.096 | K |
| Critical density | ρc | 322.0 | kg/m³ |
| Critical pressure | pc | 22.064 | MPa |
| Specific gas constant | R | 0.46151805 | kJ/(kg·K) |

---

## 2. Theoretical Foundation

### 2.1 Helmholtz Free Energy Formulation

The core of IAPWS-95 is the dimensionless Helmholtz free energy:

```
φ(δ, τ) = φ°(δ, τ) + φʳ(δ, τ)
```

Where:
- **φ°** — Ideal gas part (depends only on temperature)
- **φʳ** — Residual part (accounts for intermolecular interactions)

### Reduced Variables

```
δ = ρ / ρc          (reduced density)
τ = Tc / T          (inverse reduced temperature)
```

### 2.2 Ideal Gas Part φ°(δ, τ)

Based on IAPWS-95 Equation 5 and Tables 1, 4:

```
φ°(δ, τ) = ln(δ) + n₁ + n₂τ + n₃ln(τ) + Σᵢ₌₄⁸ nᵢ·ln[1 − exp(−γᵢτ)]
```

- 8 coefficients `nᵢ` (i=1..8) from Table 1
- 5 coefficients `γᵢ` (i=4..8) from Table 4
- 5 exponential terms representing quantum energy levels

### 2.3 Residual Part φʳ(δ, τ)

Based on IAPWS-95 Equation 6 and Table 5. The residual part consists of **56 terms** in four categories:

| Category | Index Range | Count | Form |
|----------|-------------|-------|------|
| Polynomial | i = 1–7 | 7 | nᵢ·δᵈⁱ·τᵗⁱ |
| Exponential (c=1) | i = 8–22 | 15 | nᵢ·δᵈⁱ·τᵗⁱ·exp(−δ) |
| Exponential (c=2) | i = 23–42 | 20 | nᵢ·δᵈⁱ·τᵗⁱ·exp(−δ²) |
| Exponential (c=3,4,6) | i = 43–51 | 9 | nᵢ·δᵈⁱ·τᵗⁱ·exp(−δᶜⁱ) |
| Gaussian | i = 52–54 | 3 | nᵢ·δᵈⁱ·τᵗⁱ·exp[−αᵢ(δ−εᵢ)² − βᵢ(τ−γᵢ)²] |
| Non-analytic | i = 55–56 | 2 | nᵢ·Δᵇⁱ·δF(δ, τ) |

The non-analytic terms model the critical region near the critical point using a non-analytic function F(δ, τ).

### 2.4 Thermodynamic Property Formulas

All properties are derived from φ and its partial derivatives (IAPWS-95 Table 3):

| Property | Formula | Unit |
|----------|---------|------|
| Pressure p | R·T·ρ·(1 + δ·∂φʳ/∂δ) / 1000 | MPa |
| Internal energy u | R·T·τ·(∂φ°/∂τ + ∂φʳ/∂τ) | kJ/kg |
| Enthalpy h | R·T·[τ·(∂φ°/∂τ + ∂φʳ/∂τ) + 1 + δ·∂φʳ/∂δ] | kJ/kg |
| Entropy s | R·[τ·(∂φ°/∂τ + ∂φʳ/∂τ) − φ° − φʳ] | kJ/(kg·K) |
| Isochoric heat capacity cv | −R·τ²·(∂²φ°/∂τ² + ∂²φʳ/∂τ²) | kJ/(kg·K) |
| Isobaric heat capacity cp | cv + R·N²/D | kJ/(kg·K) |
| Speed of sound w | √[R·T·(1 + 2δ·∂φʳ/∂δ + δ²·∂²φʳ/∂δ² − N²/(τ²·∂²φ/∂τ²))] · √1000 | m/s |

Where:
```
N = 1 + δ·∂φʳ/∂δ − δ·τ·∂²φʳ/∂δ∂τ
D = 1 + 2δ·∂φʳ/∂δ + δ²·∂²φʳ/∂δ²
```

---

## 3. Project Structure

```
iapws95_rust/
├── Cargo.toml                          # Package manifest (no external deps)
├── src/
│   ├── lib.rs                          # Library entry point, re-exports modules
│   ├── iapws95.rs                      # Main module: constants, API functions
│   ├── iapws95_ideal.rs                # Ideal gas part φ° and derivatives
│   ├── iapws95_residual.rs             # Residual part φʳ and derivatives
│   └── iapws95_saturation.rs           # Saturation properties solver
├── examples/
│   └── basic_usage.rs                  # Usage example (single-phase + saturation)
└── tests/
    ├── td_free_energy.rs               # Helmholtz free energy verification (Table 6)
    ├── td_test.rs                      # T-ρ-p EOS verification (Table 7)
    └── T_saturation_table8.rs          # Saturation properties verification (Table 8)
```

### Module Dependencies

```
lib.rs
├── iapws95.rs              → depends on: ideal, residual
├── iapws95_ideal.rs        → no dependencies
├── iapws95_residual.rs     → no dependencies
└── iapws95_saturation.rs   → depends on: iapws95 (which includes ideal + residual)
```

---

## 4. Module Reference

### 4.1 `iapws95` — Main Module

**File:** `src/iapws95.rs`

This module provides reference constants, helper functions for reduced variables, and the main thermodynamic property calculation API.

#### Constants

| Name | Value | Description |
|------|-------|-------------|
| `IAPWS95_TCRIT` | 647.096 | Critical temperature (K) |
| `IAPWS95_RHOCRIT` | 322.0 | Critical density (kg/m³) |
| `IAPWS95_PCRIT` | 22.064 | Critical pressure (MPa) |
| `IAPWS95_R` | 0.46151805 | Specific gas constant (kJ/(kg·K)) |
| `IAPWS95_TMIN` | 273.16 | Minimum temperature (K) |
| `IAPWS95_TMAX` | 1273.0 | Maximum temperature (K) |
| `IAPWS95_PMAX` | 1000.0 | Maximum pressure (MPa) |

#### Helper Functions

##### `reduced_density(rho: f64) -> f64`

Calculates the reduced density δ = ρ/ρc.

```rust
pub fn reduced_density(rho: f64) -> f64 {
    rho / IAPWS95_RHOCRIT
}
```

**Parameters:**
- `rho` — Density in kg/m³

**Returns:** Reduced density δ (dimensionless)

##### `inv_reduced_temp(T: f64) -> f64`

Calculates the inverse reduced temperature τ = Tc/T.

```rust
pub fn inv_reduced_temp(T: f64) -> f64 {
    IAPWS95_TCRIT / T
}
```

**Parameters:**
- `T` — Temperature in K

**Returns:** Inverse reduced temperature τ (dimensionless)

#### Main API Functions (Kelvin input)

All functions accept temperature in Kelvin and density in kg/m³.

##### `calc_pressure(T: f64, rho: f64) -> f64`

Calculates pressure from the equation of state:

```
p = R·T·ρ·(1 + δ·∂φʳ/∂δ) / 1000    [MPa]
```

##### `calc_internal_energy(T: f64, rho: f64) -> f64`

Calculates specific internal energy:

```
u = R·T·τ·(∂φ°/∂τ + ∂φʳ/∂τ)    [kJ/kg]
```

##### `calc_enthalpy(T: f64, rho: f64) -> f64`

Calculates specific enthalpy:

```
h = R·T·[τ·(∂φ°/∂τ + ∂φʳ/∂τ) + 1 + δ·∂φʳ/∂δ]    [kJ/kg]
```

##### `calc_entropy(T: f64, rho: f64) -> f64`

Calculates specific entropy:

```
s = R·[τ·(∂φ°/∂τ + ∂φʳ/∂τ) − φ° − φʳ]    [kJ/(kg·K)]
```

##### `calc_cv(T: f64, rho: f64) -> f64`

Calculates constant-volume specific heat:

```
cv = −R·τ²·(∂²φ°/∂τ² + ∂²φʳ/∂τ²)    [kJ/(kg·K)]
```

##### `calc_cp(T: f64, rho: f64) -> f64`

Calculates constant-pressure specific heat:

```
cp = cv + R·(1 + δ·∂φʳ/∂δ − δ·τ·∂²φʳ/∂δ∂τ)² / (1 + 2δ·∂φʳ/∂δ + δ²·∂²φʳ/∂δ²)    [kJ/(kg·K)]
```

##### `calc_speed_of_sound(T: f64, rho: f64) -> f64`

Calculates speed of sound:

```
w = √[R·T·(1 + 2δ·∂φʳ/∂δ + δ²·∂²φʳ/∂δ² − N²/(τ²·∂²φ/∂τ²))] · √1000    [m/s]
```

Where `N = 1 + δ·∂φʳ/∂δ − δ·τ·∂²φʳ/∂δ∂τ`. The factor √1000 converts from kJ/kg to J/kg.

#### Convenience Functions (Celsius input)

These functions convert Celsius to Kelvin and delegate to the Kelvin-based API:

| Function | Description |
|----------|-------------|
| `tr2p(t_c, rho)` | Pressure at T(°C), ρ |
| `tr2u(t_c, rho)` | Internal energy at T(°C), ρ |
| `tr2h(t_c, rho)` | Enthalpy at T(°C), ρ |
| `tr2s(t_c, rho)` | Entropy at T(°C), ρ |
| `tr2cv(t_c, rho)` | cv at T(°C), ρ |
| `tr2cp(t_c, rho)` | cp at T(°C), ρ |
| `tr2w(t_c, rho)` | Speed of sound at T(°C), ρ |

Implementation pattern:
```rust
pub fn tr2p(t_c: f64, rho: f64) -> f64 {
    let t_k = t_c + 273.15;
    calc_pressure(t_k, rho)
}
```

#### Range Check

##### `iapws95_in_range(T: f64, _p: Option<f64>) -> bool`

Checks if the temperature is within the valid range [273.16 K, 1273 K].

---

### 4.2 `iapws95_ideal` — Ideal Gas Part

**File:** `src/iapws95_ideal.rs`

Implements φ°(δ, τ) and all its partial derivatives up to second order.

#### Coefficients

| Array | Size | Description |
|-------|------|-------------|
| `IDEAL_N` | 8 | n₁ through n₈ (Table 1 of IAPWS-95) |
| `IDEAL_GAMMA` | 5 | γ₄ through γ₈ (Table 4 of IAPWS-95) |

#### Functions

##### `phi_ideal(delta: f64, tau: f64) -> f64`

Computes φ°(δ, τ). The implementation follows the formula directly with a loop over the 5 exponential terms.

##### `dphi_ideal_ddelta(delta: f64) -> f64`

Returns 1/δ (analytical first derivative w.r.t. δ).

##### `d2phi_ideal_ddelta2(delta: f64) -> f64`

Returns −1/δ² (analytical second derivative w.r.t. δ).

##### `dphi_ideal_dtau(tau: f64) -> f64`

Computes ∂φ°/∂τ using the analytical derivative formula with a loop over exponential terms.

##### `d2phi_ideal_dtau2(tau: f64) -> f64`

Computes ∂²φ°/∂τ² using the analytical second derivative formula.

##### `d2phi_ideal_dtaudelta(_delta: f64, _tau: f64) -> f64`

Returns 0.0 (the ideal gas part has no δ-τ coupling).

---

### 4.3 `iapws95_residual` — Residual Part

**File:** `src/iapws95_residual.rs`

Implements φʳ(δ, τ) and all its partial derivatives up to second order for the 56-term formulation.

#### Data Structures

| Struct | Fields | Purpose |
|--------|--------|---------|
| `PolyTerm` | d, t, n | Polynomial term: δᵈ·τᵗ·n |
| `ExpTermC1C2` | d, t, n | Exponential term with c=1 or c=2 |
| `ExpTermCN` | c, d, t, n | Exponential term with variable c (3, 4, or 6) |
| `GaussTerm` | d, t, n, a, b, g, e | Gaussian term: exp[−α(δ−ε)² − β(τ−γ)²] |
| `NonAnalTerm` | a, b, B, n, C, D, A, bt | Non-analytic critical region term |

#### Coefficient Arrays

| Array | Size | Category | Index Mapping |
|-------|------|----------|---------------|
| `RES_POLY_D1` | 7 | Polynomial (i=1–7) | Direct mapping |
| `RES_EXP_D2_C1` | 15 | Exponential c=1 (i=8–22) | i = idx + 1 |
| `RES_EXP_D2_C2` | 20 | Exponential c=2 (i=23–42) | i = idx + 1 |
| `RES_EXP_D2_CN` | 9 | Exponential c=3,4,6 (i=43–51) | Direct mapping |
| `RES_GAUSS` | 3 | Gaussian (i=52–54) | Direct mapping |
| `RES_NON_ANAL` | 2 | Non-analytic (i=55–56) | Direct mapping |

#### Functions

##### `phi_residual(delta: f64, tau: f64) -> f64`

Computes the full residual Helmholtz free energy by summing all 56 terms across all categories. Precomputes delta powers (δ², δ³, δ⁴, δ⁶) for efficiency in exponential term evaluation.

##### `dphi_residual_ddelta(delta: f64, tau: f64) -> f64`

Computes ∂φʳ/∂δ using analytical derivatives of each term category. For non-analytic terms, applies the chain rule to compute both the direct δ derivative and the indirect contribution through Δ(δ, τ).

##### `d2phi_residual_ddelta2(delta: f64, tau: f64) -> f64`

Computes ∂²φʳ/∂δ². Polynomial terms with d < 2 contribute zero to the second derivative.

##### `dphi_residual_dtau(delta: f64, tau: f64) -> f64`

Computes ∂φʳ/∂τ using analytical derivatives.

##### `d2phi_residual_dtau2(delta: f64, tau: f64) -> f64`

Computes ∂²φʳ/∂τ².

##### `d2phi_residual_ddelta_dtau(delta: f64, tau: f64) -> f64`

Computes the mixed derivative ∂²φʳ/∂δ∂τ.

---

### 4.4 `iapws95_saturation` — Saturation Properties

**File:** `src/iapws95_saturation.rs`

Implements saturation property calculations along the vapor-liquid equilibrium line using a hybrid approach combining IAPWS SR1-86 explicit equations with Newton's method for IAPWS-95 phase equilibrium.

#### Data Structures

##### `SaturationProperties`

| Field | Description | Unit |
|-------|-------------|------|
| `p_sat` | Saturation vapor pressure | MPa |
| `rho_l` | Saturated liquid density | kg/m³ |
| `rho_v` | Saturated vapor density | kg/m³ |
| `h_l` | Saturated liquid enthalpy | kJ/kg |
| `h_v` | Saturated vapor enthalpy | kJ/kg |
| `s_l` | Saturated liquid entropy | kJ/(kg·K) |
| `s_v` | Saturated vapor entropy | kJ/(kg·K) |

#### Algorithm Details

**Step 1: Initial Guesses (IAPWS SR1-86)**

Saturated liquid density using Wagner-type correlation:
```
ρ'/ρc = 1 + b₁τ^(1/3) + b₂τ^(2/3) + b₃τ^(5/3) + b₄τ^(16/3) + b₅τ^(43/3) + b₆τ^(110/3)
```

Saturated vapor density using exponential correlation:
```
ln(ρ''/ρc) = c₁τ^(2/6) + c₂τ^(4/6) + c₃τ^(8/6) + c₄τ^(18/6) + c₅τ^(37/6) + c₆τ^(71/6)
```

**Step 2: Phase Equilibrium Refinement (Newton's Method)**

Solves the system of two equations for δ' and δ'':

```
F1(δL, δV) = K(δV, τ) − K(δL, τ) = 0    (equal chemical potential)
F2(δL, δV) = J(δV, τ) − J(δL, τ) = 0    (equal pressure)
```

Where:
- `J(δ, τ) = δ·[1 + δ·∂φʳ/∂δ]` — dimensionless pressure term
- `K(δ, τ) = δ·∂φʳ/∂δ + φʳ + ln(δ)` — dimensionless chemical potential term

The Jacobian matrix:
```
| −∂K/∂δL    ∂K/∂δV |
| −∂J/∂δL    ∂J/∂δV |
```

**Step 3: Property Calculation**

Once δ' and δ'' are determined, all saturation properties are computed using the main IAPWS-95 API functions.

#### Functions

##### `calc_saturation_properties(T: f64) -> Option<SaturationProperties>`

Computes saturation properties at a given temperature in Kelvin. Returns `None` if T is outside [273.16 K, 647.096 K] or if the Newton solver fails to converge.

##### `sat_t(t_c: f64) -> Option<SaturationProperties>`

Same as above but accepts temperature in Celsius. Returns `None` if t_c is outside [0.01°C, 373.946°C].

---

## 5. API Reference Summary

### Single-Phase Properties (Kelvin input)

```rust
use iapws95::iapws95::*;

// Input: T in Kelvin, rho in kg/m³
let p    = calc_pressure(T, rho);        // MPa
let u    = calc_internal_energy(T, rho); // kJ/kg
let h    = calc_enthalpy(T, rho);        // kJ/kg
let s    = calc_entropy(T, rho);         // kJ/(kg·K)
let cv   = calc_cv(T, rho);              // kJ/(kg·K)
let cp   = calc_cp(T, rho);              // kJ/(kg·K)
let w    = calc_speed_of_sound(T, rho);  // m/s
```

### Single-Phase Properties (Celsius input)

```rust
use iapws95::iapws95::*;

// Input: t_c in Celsius, rho in kg/m³
let p    = tr2p(t_c, rho);        // MPa
let u    = tr2u(t_c, rho);        // kJ/kg
let h    = tr2h(t_c, rho);        // kJ/kg
let s    = tr2s(t_c, rho);        // kJ/(kg·K)
let cv   = tr2cv(t_c, rho);       // kJ/(kg·K)
let cp   = tr2cp(t_c, rho);       // kJ/(kg·K)
let w    = tr2w(t_c, rho);        // m/s
```

### Saturation Properties

```rust
use iapws95::iapws95_saturation::*;

// Kelvin input
if let Some(sat) = calc_saturation_properties(T_k) {
    println!("p_sat: {} MPa", sat.p_sat);
    println!("rho_l: {} kg/m³", sat.rho_l);
    println!("rho_v: {} kg/m³", sat.rho_v);
    println!("h_l: {} kJ/kg", sat.h_l);
    println!("h_v: {} kJ/kg", sat.h_v);
    println!("s_l: {} kJ/(kg·K)", sat.s_l);
    println!("s_v: {} kJ/(kg·K)", sat.s_v);
}

// Celsius input
if let Some(sat) = sat_t(t_celsius) { ... }
```

### Reduced Variables

```rust
use iapws95::iapws95::*;

let delta = reduced_density(rho);   // δ = ρ/ρc
let tau   = inv_reduced_temp(T);    // τ = Tc/T
```

---

## 6. Saturation Properties Algorithm — Detailed Analysis

### 6.1 Why a Hybrid Approach?

The vapor-liquid equilibrium problem is inherently nonlinear and can be sensitive to initial guesses, especially near the critical point. The hybrid approach leverages:

1. **SR1-86 explicit equations** — Provide robust initial estimates across the entire valid temperature range with typical errors < 0.5% for densities
2. **Newton's method on IAPWS-95** — Refines the solution to machine precision using the more accurate IAPWS-95 formulation

### 6.2 Convergence Strategy

The Newton solver includes several safeguards:

| Safeguard | Description |
|-----------|-------------|
| Density clamping | δL ∈ [1.0, 4.5], δV ∈ [1e−7, 0.8] |
| Damping | Halved iteratively if new densities violate bounds |
| Bounds check | δL > 1.1·δV (liquid must be denser than vapor) |
| Early exit | Checks residual < 1e−12 after each iteration |
| Max iterations | 100 Newton steps + up to 20 damping attempts per step |

### 6.3 Phase Equilibrium Conditions

The solver enforces the **Maxwell equal-area criterion**:

```
p(δ', τ) = p(δ'', τ)     (mechanical equilibrium)
μ(δ', τ) = μ(δ'', τ)     (chemical equilibrium)
```

These conditions ensure that liquid and vapor phases coexist in thermodynamic equilibrium at the saturation temperature.

---

## 7. Test Suite

### 7.1 Test Files

| File | Purpose | Reference Data |
|------|---------|----------------|
| `td_free_energy.rs` | Verifies φ° and φʳ calculation accuracy | IAPWS-95 Table 6 |
| `td_test.rs` | Verifies T-ρ-p equation of state and derived properties | IAPWS-95 Table 7 (11 test points) |
| `T_saturation_table8.rs` | Verifies saturation property calculations | IAPWS-95 Table 8 |

### 7.2 Test Methodology

All tests use the `assert_approx_eq!` macro from the `assert_approx_eq` crate for floating-point comparison with specified tolerances:

| Property | Typical Tolerance |
|----------|-------------------|
| Helmholtz free energy φ° | 1e−6 |
| Helmholtz free energy φʳ | 1e−6 |
| Pressure p | 1e−4 to 1e−6 |
| Specific heat cv | 1e−4 |
| Speed of sound w | 1e−5 |
| Entropy s | 1e−6 |

### 7.3 Example Test Case (Table 7)

```rust
// IAPWS-95 Table 7: T = 300 K, ρ = 996.556 kg/m³
let T_test = 300.0;
let rho_test = 996.556;

assert_approx_eq!(0.001034827, calc_pressure(T_test, rho_test), 1e-6);
assert_approx_eq!(4.177450, calc_cv(T_test, rho_test), 1e-4);
assert_approx_eq!(1497.0, calc_speed_of_sound(T_test, rho_test), 1e-5);
```

---

## 8. Build and Usage

### 8.1 Prerequisites

- Rust toolchain (edition 2021)
- Cargo build system

### 8.2 Building

```bash
# Debug build
cargo build

# Release build (optimized)
cargo build --release

# Run tests
cargo test

# Generate and open documentation
cargo doc --open
```

### 8.3 Running Examples

```bash
cargo run --example basic_usage
```

### 8.4 Adding as a Dependency

Add to your `Cargo.toml`:

```toml
[dependencies]
iapws95 = { git = "https://github.com/thermalogic/iapws95_rust" }
```

Or if published to crates.io:

```toml
[dependencies]
iapws95 = "0.1"
```

### 8.5 Adding New Properties

To add a new thermodynamic property based on Helmholtz free energy derivatives:

```rust
// In src/iapws95.rs
pub fn calc_new_property(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    
    // Access derivatives from iapws95_ideal and iapws95_residual
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    // ... apply the appropriate formula
    
    result
}
```

---

## Appendix A: Term Classification Summary

### Residual Part — 56 Terms Breakdown

| Category | Count | c Value | δ Power Range | τ Power Range |
|----------|-------|---------|---------------|---------------|
| Polynomial | 7 | — | 1–4 | −0.5 to 1.0 |
| Exponential (c=1) | 15 | 1 | 1–15 | 1–13 |
| Exponential (c=2) | 20 | 2 | 1–12 | 1–10 |
| Exponential (c=3) | 4 | 3 | 3–5 | 16–23 |
| Exponential (c=4) | 1 | 4 | 14 | 10 |
| Exponential (c=6) | 4 | 6 | 3–6 | 44–50 |
| Gaussian | 3 | — | 3 | 0, 1, 4 |
| Non-analytic | 2 | — | — | — |

## Appendix B: Unit Conventions

| Quantity | SI Base | Library Output | Notes |
|----------|---------|----------------|-------|
| Temperature | K | K (input) or °C (convenience API) | Kelvin is the canonical unit |
| Density | kg/m³ | kg/m³ | Direct input/output |
| Pressure | Pa | MPa | Divided by 1e6 in output |
| Energy | J/kg | kJ/kg | Divided by 1e3 in output |
| Entropy | J/(kg·K) | kJ/(kg·K) | Divided by 1e3 in output |
| Speed of sound | m/s | m/s | No conversion needed |

---

*This documentation corresponds to IAPWS-95 Formulation 1995 (Revised 2018). For the official standard, refer to [IAPWS R6-95(2018)](https://iapws.org/readme/iapws-r1/).*
