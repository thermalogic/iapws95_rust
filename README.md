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
let t_c = 100.0;      // Temperature in Celsius
let rho = 958.0;       // Density (kg/m³) — saturated water at 100°C

let p    = tr2p(t_c, rho);        // Pressure (MPa)
let h    = tr2h(t_c, rho);        // Enthalpy (kJ/kg)
let s    = tr2s(t_c, rho);        // Entropy (kJ/(kg·K))
let w    = tr2w(t_c, rho);        // Speed of sound (m/s)
```

## API Reference

### Main Functions (Temperature in °C)

These convenience functions accept temperature in **Celsius** and provide shorter function names as the primary API:

| Function | Description | Returns |
|------|------|------|
| `tr2p(t_c, rho)` | Pressure at T(°C), ρ(kg/m³) | MPa |
| `tr2u(t_c, rho)` | Internal energy at T(°C), ρ(kg/m³) | kJ/kg |
| `tr2h(t_c, rho)` | Enthalpy at T(°C), ρ(kg/m³) | kJ/kg |
| `tr2s(t_c, rho)` | Entropy at T(°C), ρ(kg/m³) | kJ/(kg·K) |
| `tr2cv(t_c, rho)` | Constant-volume specific heat at T(°C), ρ(kg/m³) | kJ/(kg·K) |
| `tr2cp(t_c, rho)` | Constant-pressure specific heat at T(°C), ρ(kg/m³) | kJ/(kg·K) |
| `tr2w(t_c, rho)` | Speed of sound at T(°C), ρ(kg/m³) | m/s |

### Main Functions (Temperature in K)

These functions accept temperature in **Kelvin**:

| Function | Description | Returns |
|------|------|------|
| `calc_pressure(T, rho)` | Pressure at T(K), ρ(kg/m³) | MPa |
| `calc_internal_energy(T, rho)` | Internal energy at T(K), ρ(kg/m³) | kJ/kg |
| `calc_enthalpy(T, rho)` | Enthalpy at T(K), ρ(kg/m³) | kJ/kg |
| `calc_entropy(T, rho)` | Entropy at T(K), ρ(kg/m³) | kJ/(kg·K) |
| `calc_cv(T, rho)` | Constant-volume specific heat at T(K), ρ(kg/m³) | kJ/(kg·K) |
| `calc_cp(T, rho)` | Constant-pressure specific heat at T(K), ρ(kg/m³) | kJ/(kg·K) |
| `calc_speed_of_sound(T, rho)` | Speed of sound at T(K), ρ(kg/m³) | m/s |

### Saturation Properties

```rust
use iapws95::iapws95_saturation::*;

if let Some(props) = sat_t(100.0) {  // Temperature in °C
    println!("p_sat: {} MPa", props.p_sat);
    println!("rho_l: {} kg/m³", props.rho_l);
    println!("rho_v: {} kg/m³", props.rho_v);
}
```

## Testing

```bash
cargo test                    # Run all tests
cargo test --test td_test    # Specific test
```

| Test | Description |
|------|------|
| `td_free_energy.rs` | Helmholtz free energy verification (Table 6) |
| `td_test.rs` | T-ρ-p equation of state (Table 7) |
| `T_saturation_table8.rs` | Saturation properties (Table 8) |

## Dependencies

No external runtime dependencies. Only uses Rust standard library.

```toml
[dev-dependencies]
assert_approx_eq = "1.1.0"
```

## References

1. [IAPWS R6-95(2018)](https://iapws.org/documents/release/IAPWS-95/) - Revised Release on the IAPWS Formulation 1995
