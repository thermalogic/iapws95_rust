//! IAPWS-95 Thermodynamic Properties Library for Water/Steam
//!
//! Rust implementation of the International Association for the Properties
//! of Water and Steam Formulation 1995 (Revised 2018)
//!
//! Range: 273.16 K to 1273 K, up to 1000 MPa
//!
//! # Features
//!
//! - `python` - Enable Python bindings via PyO3
//! - `cffi` - Enable C FFI bindings (default for non-test builds)
//! - `all-bindings` - Enable all optional bindings

pub mod iapws95;
pub mod iapws95_ideal;
pub mod iapws95_residual;
pub mod iapws95_saturation;

#[cfg(feature = "python")]
pub mod py_iapws95;

#[cfg(any(feature = "cffi", not(test)))]
pub mod c_iapws95;
