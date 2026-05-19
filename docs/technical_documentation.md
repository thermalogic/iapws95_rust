# IAPWS-95 Thermodynamic Properties Library — Technical Documentation

**Language:** Rust (Edition 2021)  
**Standard:** IAPWS-95 (Revised 2018)
**Version:** 0.2.1

## Table of Contents

1. [Overview](#1-overview)
   - [Installation](#installation)
   - [Key Features](#key-features)
   - [Valid Range](#valid-range)
   - [Reference Constants](#reference-constants)
2. [Theoretical Foundation](#2-theoretical-foundation)
3. [Project Structure](#3-project-structure)
4. [Module Reference](#4-module-reference)
5. [API Reference](#5-api-reference)
6. [Saturation Properties Algorithm](#6-saturation-properties-algorithm)
7. [Test Suite](#7-test-suite)
8. [Build and Usage](#8-build-and-usage)
9. [Performance Optimizations](#9-performance-optimizations)
10. [References](#references)

## 1. Overview

This library implements the **IAPWS-95** formulation — the International Association for the Properties of Water and Steam standard for thermodynamic properties of ordinary water substance. It provides accurate calculations of pressure, enthalpy, entropy, specific heats, speed of sound, and saturation properties for water and steam.

### Installation

The package is published on [crates.io](https://crates.io/crates/iapws95). Install it via:

```bash
cargo install iapws95
```

Or add to your `Cargo.toml`:

```toml
[dependencies]
iapws95 = "0.2"
```

### Key Features

- Zero external runtime dependencies (Rust standard library only)
- Full IAPWS-95 Helmholtz free energy formulation (56 terms across 6 categories)
- Single-phase property calculations in the (T, ρ) domain
- Saturation properties with phase equilibrium solver
- Verified against IAPWS reference tables (Tables 6, 7, 8)
- Optimized with LTO and compiler optimizations for production use

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

Based on IAPWS-95 Equation 6 and Table 5. The residual part consists of **56 terms** in six categories:

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
| Joule-Thomson coefficient μ | δ(τ·φʳ_δτ − φʳ_δ) / [ρ·Cp·(1 + 2δφʳ_δ + δ²φʳ_δδ)] | K/MPa |
| Isothermal throttling coefficient (∂τ/∂p)_T | 1 − (1 + δφʳ_δ − δτφʳ_δτ) / (1 + 2δφʳ_δ + δ²φʳ_δδ) | kJ/(kg·MPa) |
| Isentropic temperature-pressure coefficient β_s | (1 + δφʳ_δ − δτφʳ_δτ) / [ρR((1 + δφʳ_δ − δτφʳ_δτ)² − τ²(φ°_ττ + φʳ_ττ)(1 + 2δφʳ_δ + δ²φʳ_δδ))] | 1/K |

Where:
```
N = 1 + δ·∂φʳ/∂δ − δ·τ·∂²φʳ/∂δ∂τ
D = 1 + 2δ·∂φʳ/∂δ + δ²·∂²φʳ/∂δ²
```

## 3. Project Structure

```
iapws95_rust/
├── Cargo.toml                          # Package manifest (no external deps)
├── .cargo/
│   └── config.toml                     # Build configuration (compiler flags)
├── src/
│   ├── lib.rs                          # Library entry point, re-exports modules
│   ├── iapws95.rs                      # Main module: constants, API functions
│   ├── iapws95_ideal.rs                # Ideal gas part φ° and derivatives
│   ├── iapws95_residual.rs             # Residual part φʳ and derivatives
│   └── iapws95_saturation.rs           # Saturation properties solver
├── benches/
│   └── iapws95_bench.rs                # Performance benchmarks (Criterion.rs)
├── examples/
│   └── basic_usage.rs                  # Usage example (single-phase + saturation)
└── tests/
    ├── td_free_energy.rs               # Helmholtz free energy verification (Table 6)
    ├── td_test_table7.rs               # T-ρ-p EOS verification (Table 7)
    ├── td_test_95_97.rs                # IAPWS-95 vs IAPWS-IF97 comparison
    ├── T_saturation_table8.rs          # Saturation properties verification (Table 8)
    └── td_data.rs                      # Shared test data module
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

#### Public API Functions (Celsius input)

These are the primary public API functions. They accept temperature in **°C** and density in kg/m³:

| Function          | Description                                    | Returns    |
| ----------------- | ---------------------------------------------- | ---------- |
| `tr2p(t_c, rho)`  | Pressure                                       | MPa        |
| `tr2u(t_c, rho)`  | Internal energy                                | kJ/kg      |
| `tr2h(t_c, rho)`  | Enthalpy                                       | kJ/kg      |
| `tr2s(t_c, rho)`  | Entropy                                        | kJ/(kg·K)  |
| `tr2cv(t_c, rho)` | Constant-volume specific heat                  | kJ/(kg·K)  |
| `tr2cp(t_c, rho)` | Constant-pressure specific heat                | kJ/(kg·K)  |
| `tr2w(t_c, rho)`  | Speed of sound                                 | m/s        |
| `tr2jt(t_c, rho)` | Joule-Thomson coefficient                      | K/MPa      |
| `tr2itt(t_c, rho)`| Isothermal throttling coefficient              | kJ/(kg·MPa)|
| `tr2beta_s(t_c, rho)` | Isentropic temperature-pressure coefficient| 1/K        |


Implementation pattern:
```rust
pub fn tr2p(t_c: f64, rho: f64) -> f64 {
    let t_k = t_c + 273.15;
    calc_pressure(t_k, rho)
}
```

#### Internal Functions (Kelvin input)

These functions are `pub(crate)` and used internally by the `tr2*` convenience functions:

| Function | Description |
|----------|-------------|
| `calc_pressure(T, rho)` | Pressure at T(K), ρ |
| `calc_internal_energy(T, rho)` | Internal energy at T(K), ρ |
| `calc_enthalpy(T, rho)` | Enthalpy at T(K), ρ |
| `calc_entropy(T, rho)` | Entropy at T(K), ρ |
| `calc_cv(T, rho)` | cv at T(K), ρ |
| `calc_cp(T, rho)` | cp at T(K), ρ |
| `calc_speed_of_sound(T, rho)` | Speed of sound at T(K), ρ |
| `calc_joule_thomson(T, rho)` | Joule-Thomson coefficient at T(K), ρ |
| `calc_isothermal_throttling(T, rho)` | Isothermal throttling coefficient at T(K), ρ |
| `calc_isentropic_temp_pressure(T, rho)` | Isentropic temperature-pressure coefficient at T(K), ρ |

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

Computes the full residual Helmholtz free energy by summing all 56 terms across all categories. Precomputes exponential factors (exp(−δ), exp(−δ²), exp(−δ³), exp(−δ⁴), exp(−δ⁶)) once per function call, reusing them across all terms in each category.

##### `dphi_residual_ddelta(delta: f64, tau: f64) -> f64`

Computes ∂φʳ/∂δ using analytical derivatives of each term category. Precomputes exponential factors for efficiency. For non-analytic terms, applies the chain rule to compute both the direct δ derivative and the indirect contribution through Δ(δ, τ).

##### `d2phi_residual_ddelta2(delta: f64, tau: f64) -> f64`

Computes ∂²φʳ/∂δ² using analytical second derivatives for each term category. Precomputes exponential factors. For terms where d < 2, uses `powi` for negative exponents:

- **Polynomial terms**: nᵢ·dᵢ·(dᵢ−1)·δ^(dᵢ−2)·τ^tᵢ (terms with d < 2 contribute zero)
- **Exponential terms (c=1)**: nᵢ·exp(−δ)·δ^(d−2)·τ^t·[(d−δ)(d−1−δ)−δ]
- **Exponential terms (c=2)**: nᵢ·exp(−δ²)·δ^(d−2)·τ^t·[(d−1−2δ²)(d−2δ²)−4δ²]
- **Exponential terms (c=3,4,6)**: nᵢ·exp(−δᶜ)·δ^(d−2)·τ^t·[(d−1−c·δᶜ)(d−c·δᶜ)−c²·δᶜ]
- **Gaussian terms**: Uses product rule with first and second derivative factors
- **Non-analytic terms**: Full chain rule implementation with Δ(δ,τ) second derivatives

##### `dphi_residual_dtau(delta: f64, tau: f64) -> f64`

Computes ∂φʳ/∂τ using analytical derivatives. Precomputes exponential factors.

##### `d2phi_residual_dtau2(delta: f64, tau: f64) -> f64`

Computes ∂²φʳ/∂τ². Precomputes exponential factors.

##### `d2phi_residual_ddelta_dtau(delta: f64, tau: f64) -> f64`

Computes the mixed derivative ∂²φʳ/∂δ∂τ. Precomputes exponential factors.

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
let mu   = calc_joule_thomson(T, rho);   // K/MPa
let itt  = calc_isothermal_throttling(T, rho); // dimensionless
let beta = calc_isentropic_temp_pressure(T, rho); // 1/K
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
let mu   = tr2jt(t_c, rho);       // K/MPa
let itt  = tr2itt(t_c, rho);      // dimensionless
let beta = tr2beta_s(t_c, rho);   // 1/K
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
| `td_test_table7.rs` | Verifies T-ρ-p equation of state and derived properties | IAPWS-95 Table 7 (11 test points) |
| `td_test_95_97.rs` | Compares IAPWS-95 vs IAPWS-IF97 calculations | seuif97 library |
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

### Building

```bash
cargo build --release
```

The release profile is optimized with:
- `opt-level = 3` - Maximum optimization level
- `lto = "thin"` - Thin Link Time Optimization for cross-module optimization
- `codegen-units = 1` - Single codegen unit for better optimization

### Usage

```rust
use iapws95::iapws95::*;

// Single-phase properties (Celsius input)
let t_c = 26.85;  // 300 K
let rho = 996.556;

let p = tr2p(t_c, rho);  // Pressure in MPa
let h = tr2h(t_c, rho);  // Enthalpy in kJ/kg
let s = tr2s(t_c, rho);  // Entropy in kJ/(kg·K)

// Saturation properties
use iapws95::iapws95_saturation::sat_t;

if let Some(sat) = sat_t(100.0) {
    println!("Saturation pressure: {} MPa", sat.p_sat);
}
```

### Testing

```bash
cargo test                    # Run all tests
cargo test --test td_test_table7  # Run specific test
```

### Benchmarking

```bash
cargo bench
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

## 9. Performance Optimizations

### 9.1 Compiler Optimizations

The release profile is configured for maximum performance:

```toml
[profile.release]
opt-level = 3           # Maximum optimization level
lto = "thin"            # Thin Link Time Optimization
codegen-units = 1       # Single codegen unit for better optimization
```

### 9.2 Algorithm Optimizations

- **Precomputed exponential factors**: Exponential terms (exp(−δ), exp(−δ²), exp(−δ³), exp(−δ⁴), exp(−δ⁶)) are computed once per function call and reused across all terms in each category
- **Inline functions**: All property calculation functions are marked with `#[inline]` for zero-overhead function calls
- **Analytical derivatives**: All derivatives are computed analytically rather than numerically, avoiding iterative approximation

### 9.3 Memory Efficiency

- **Zero heap allocation**: All calculations use stack-allocated arrays and primitive types
- **No dynamic memory**: Coefficient arrays are compile-time constants
- **Minimal stack usage**: Typical property calculation uses < 1KB stack space

### 9.4 Benchmarking

Performance benchmarks are available in `benches/iapws95_bench.rs` using Criterion.rs:

```bash
cargo bench
```

Typical performance characteristics:
- Single-phase property calculation: ~100-500 nanoseconds
- Saturation property calculation: ~5-20 microseconds (includes Newton solver iterations)

### 9.5 Exponential Factor Reuse

Exponential factors are computed once per function call and reused across all terms in each category:

| Factor | Computed Once | Used In |
|--------|---------------|---------|
| `exp(−δ)` | c=1 exponential terms (15 terms) | `RES_EXP_D2_C1` |
| `exp(−δ²)` | c=2 exponential terms (20 terms) | `RES_EXP_D2_C2` |
| `exp(−δ³)` | c=3 exponential terms (4 terms) | `RES_EXP_D2_CN[0..4]` |
| `exp(−δ⁴)` | c=4 exponential term (1 term) | `RES_EXP_D2_CN[4]` |
| `exp(−δ⁶)` | c=6 exponential terms (4 terms) | `RES_EXP_D2_CN[5..9]` |

### 9.6 Direct cp Calculation

The `calc_cp` function computes the isobaric heat capacity directly using the full formula rather than calling `calc_cv`. This avoids redundant evaluation of `d2phi_ideal_dtau2` and `d2phi_residual_dtau2`, which are expensive second-order derivative computations.

**Before (redundant):**
```rust
let cv = calc_cv(T, rho); // computes d2phi_dtau2 internally
// then computes additional cp terms
```

**After (direct):**
```rust
let cv_part = -tau * tau * (phi_o_tt + phi_r_tt); // reuse already-computed derivatives
let numerator = (1.0 + delta * dphi_ddelta - delta * tau * d2phi_ddelta_dtau).powi(2);
let denominator = 1.0 + 2.0 * delta * dphi_ddelta + delta * delta * d2phi_ddelta2;
IAPWS95_R * (cv_part + numerator / denominator)
```

### 9.7 Negative Exponent Handling

In `d2phi_residual_ddelta2`, terms with d < 2 require negative exponents (δ^(d−2)). The implementation uses `powi` directly for these cases:

```rust
let delta_d_minus_2 = delta.powi(term.d - 2);
```

### 9.8 Performance Benchmarking

The library includes comprehensive performance benchmarks using Criterion.rs framework to measure and track the performance of all thermodynamic property calculations.

#### Benchmark Structure

The benchmark suite ([benches/iapws95_bench.rs](../benches/iapws95_bench.rs)) consists of two test groups:

1. **Properties** - Tests all properties across liquid, steam, and supercritical states
2. **All Properties x 5 States** - Tests all properties across 5 different thermodynamic states

#### Running Benchmarks

```bash
cargo bench
```

Results are generated in `target/criterion/report/` as HTML reports with detailed performance analysis.

*Note: Actual performance varies with CPU architecture and compiler optimizations.*

## References

1. [IAPWS R6-95(2018)](https://iapws.org/documents/release/IAPWS-95/) - Revised Release on the IAPWS Formulation 1995
2. Wagner, W. & Pruss, A. (2002). The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance. J. Phys. Chem. Ref. Data, 31(2), 387-535.