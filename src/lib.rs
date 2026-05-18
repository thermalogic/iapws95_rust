//! IAPWS-95 Thermodynamic Properties Library for Water/Steam
//!
//! Rust implementation of the International Association for the Properties
//! of Water and Steam Formulation 1995 (Revised 2018)
//!
//! # Range of Validity
//! - Temperature: 273.16 K to 1273 K
//! - Pressure: up to 1000 MPa
//!
//! # Public API (External Users)
//!
//! ## Single-Phase Properties (from `iapws95` module)
//! All functions accept temperature in **°C** and density in **kg/m³**:
//!
//! | Function | Description | Returns |
//! |----------|-------------|---------|
//! | `tr2p(t_c, rho)` | Pressure | MPa |
//! | `tr2u(t_c, rho)` | Internal energy | kJ/kg |
//! | `tr2h(t_c, rho)` | Enthalpy | kJ/kg |
//! | `tr2s(t_c, rho)` | Entropy | kJ/(kg·K) |
//! | `tr2cv(t_c, rho)` | Constant-volume specific heat | kJ/(kg·K) |
//! | `tr2cp(t_c, rho)` | Constant-pressure specific heat | kJ/(kg·K) |
//! | `tr2w(t_c, rho)` | Speed of sound | m/s |
//! | `tr2jt(t_c, rho)` | Joule-Thomson coefficient | K/MPa |
//! | `tr2itt(t_c, rho)` | Isothermal throttling coefficient | kJ/(kg·MPa) |
//! | `tr2beta_s(t_c, rho)` | Isentropic temperature-pressure coefficient | K/MPa |
//!
//! ## Saturation Properties (from `iapws95_saturation` module)
//!
//! | Function | Description |
//! |----------|-------------|
//! | `sat_t(t_c)` | Saturation properties at T(°C) |
//!
//! Returns `Option<SaturationProperties>` containing: `p_sat`, `rho_l`, `rho_v`, `h_l`, `h_v`, `s_l`, `s_v`
//!
//! # Modules
//!
//! - `iapws95` - Main thermodynamic property calculations (10 public `tr2*` functions + internal helpers)
//! - `iapws95_ideal` - Ideal gas part of the dimensionless Helmholtz free energy φ°(δ,τ) (internal)
//! - `iapws95_residual` - Residual part of the dimensionless Helmholtz free energy φʳ(δ,τ) (internal)
//! - `iapws95_saturation` - Saturation properties calculation (1 public `sat_t` function)
//!
//! # Key Variables
//!
//! - δ (delta) = ρ/ρc - Reduced density (dimensionless)
//! - τ (tau) = Tc/T - Inverse reduced temperature (dimensionless)
//! - ρc = 322 kg/m³ - Critical density
//! - Tc = 647.096 K - Critical temperature
//!
//! # Helmholtz Free Energy Formulation
//!
//! The dimensionless Helmholtz free energy is split into two parts:
//!
//! ```text
//! φ(δ,τ) = φ°(δ,τ) + φʳ(δ,τ)
//! ```
//!
//! where φ° is the ideal gas part and φʳ is the residual part.

pub mod iapws95;
pub mod iapws95_ideal;
pub mod iapws95_residual;
pub mod iapws95_saturation;
