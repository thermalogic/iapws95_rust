# IAPWS-95 Thermodynamic Properties Library (Rust Implementation)

A Rust implementation of the [IAPWS-95](https://iapws.org/documents/release/IAPWS-95/) standard for calculating thermodynamic properties of water and steam.

## Project Structure

```
iapws95_rust/
├── src/
│   ├── lib.rs                # Library entry point
│   ├── iapws95.rs            # Main module: constants & API functions
│   ├── iapws95_ideal.rs      # Ideal gas part (φ°)
│   ├── iapws95_residual.rs   # Residual part (φʳ)
│   └── iapws95_saturation.rs # Saturation properties
├── examples/
│   └── basic_usage.rs
└── tests/
    ├── td_free_energy.rs     # Helmholtz free energy verification
    ├── td_test.rs            # T-ρ-p equation of state test
    └── T_saturation_table8.rs # Saturation properties test
```

## Quick Start

### Installation

The package is published on [crates.io](https://crates.io/crates/iapws95). Install it via:

```bash
cargo install iapws95
```

Or add to your `Cargo.toml`:

```toml
[dependencies]
iapws95 = "0.1"
```

### Build from Source

```bash
git clone https://github.com/thermalogic/iapws95_rust.git
cd iapws95_rust
cargo build --release
cargo test
cargo doc --open
```

### Usage

```rust
use iapws95::iapws95::*;
let t_c = 26.85;      // Temperature in Celsius，300.0K
let rho =  0.9965560e3;       // Density (kg/m³) 

let p    = tr2p(t_c, rho);        // Pressure (MPa)
let u    = tr2u(t_c, rho);        // Internal energy (kJ/kg)
let h    = tr2h(t_c, rho);        // Enthalpy (kJ/kg)
let s    = tr2s(t_c, rho);        // Entropy (kJ/(kg·K))
let cv   = tr2cv(t_c, rho);       // Constant-volume specific heat (kJ/(kg·K))
let cp   = tr2cp(t_c, rho);       // Constant-pressure specific heat (kJ/(kg·K))
let w    = tr2w(t_c, rho);        // Speed of sound (m/s)
```

## API Reference

### Single-Phase Properties (Temperature in °C)

All functions accept temperature in **Celsius** and density in kg/m³:

| Function          | Description                                        | Returns   |
| ----------------- | -------------------------------------------------- | --------- |
| `tr2p(t_c, rho)`  | Pressure at T(°C), ρ(kg/m³)                        | MPa       |
| `tr2u(t_c, rho)`  | Internal energy at T(°C), ρ(kg/m³)                 | kJ/kg     |
| `tr2h(t_c, rho)`  | Enthalpy at T(°C), ρ(kg/m³)                        | kJ/kg     |
| `tr2s(t_c, rho)`  | Entropy at T(°C), ρ(kg/m³)                         | kJ/(kg·K) |
| `tr2cv(t_c, rho)` | Constant-volume specific heat at T(°C), ρ(kg/m³)   | kJ/(kg·K) |
| `tr2cp(t_c, rho)` | Constant-pressure specific heat at T(°C), ρ(kg/m³) | kJ/(kg·K) |
| `tr2w(t_c, rho)`  | Speed of sound at T(°C), ρ(kg/m³)                  | m/s       |

### Saturation Properties

```rust
use iapws95::iapws95_saturation::sat_t;

if let Some(props) = sat_t(100.0) {  // Temperature in °C
    println!("p_sat: {} MPa", props.p_sat);
    println!("rho_l: {} kg/m³", props.rho_l);
    println!("rho_v: {} kg/m³", props.rho_v);
    println!("h_l: {} kJ/kg", props.h_l);
    println!("h_v: {} kJ/kg", props.h_v);
    println!("s_l: {} kJ/(kg·K)", props.s_l);
    println!("s_v: {} kJ/(kg·K)", props.s_v);
}
```

## Testing

```bash
cargo test                    # Run all tests
cargo test --test td_test    # Specific test
```

| Test                     | Description                                  |
| ------------------------ | -------------------------------------------- |
| `td_free_energy.rs`      | Helmholtz free energy verification (Table 6) |
| `td_test.rs`             | T-ρ-p equation of state (Table 7)            |
| `T_saturation_table8.rs` | Saturation properties (Table 8)              |

## Implementation Details

### Helmholtz Free Energy Derivatives

The library implements complete first and second derivatives of the Helmholtz free energy:

- **Ideal gas part (φ°)**: Analytical derivatives for all terms
- **Residual part (φʳ)**: Full derivative support including:
  - Polynomial terms (i=1-7)
  - Exponential terms with c=1,2,3,4,6 (i=8-51)
  - Gaussian terms (i=52-54)
  - Non-analytic terms (i=55-56)

All derivative formulas have been verified against reference values from IAPWS-95 Table 6.

## Dependencies

No external runtime dependencies. Only uses Rust standard library.

```toml
[dev-dependencies]
assert_approx_eq = "1.1.0"
```

## References

1. [IAPWS R6-95(2018)](https://iapws.org/documents/release/IAPWS-95/) - Revised Release on the IAPWS Formulation 1995
2. Wagner, W. & Pruss, A. (2002). The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance. J. Phys. Chem. Ref. Data, 31(2), 387-535.

