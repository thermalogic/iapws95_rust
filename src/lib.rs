//! IAPWS-95 Thermodynamic Properties Library for Water/Steam
//!
//! Rust implementation of the International Association for the Properties
//! of Water and Steam Formulation 1995 (Revised 2018)
//!
//! # Range of Validity
//! - Temperature: 273.16 K to 1273 K
//! - Pressure: up to 1000 MPa
//!
//! # Modules
//!
//! - `iapws95` - Main thermodynamic property calculations (pressure, enthalpy, entropy, speed of sound, etc.)
//! - `iapws95_ideal` - Ideal gas part of the dimensionless Helmholtz free energy φ°(δ,τ)
//! - `iapws95_residual` - Residual part of the dimensionless Helmholtz free energy φʳ(δ,τ)
//! - `iapws95_saturation` - Saturation properties calculation using hybrid SR1-86 + Newton's method
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
