//! IAPWS-95 Header - Reference Constants, Ranges, and Data Structures
//!
//! Translated from iapws95.h
//!
//! # Public API
//! - `tr2p`, `tr2u`, `tr2h`, `tr2s`, `tr2cv`, `tr2cp`, `tr2w` -`tr2jt`- `tr2itt``- `tr2beta_s`- Single-phase properties (T in °C, ρ in kg/m³)
//!
//! # Internal Functions (pub(crate))
//! - `calc_pressure`, `calc_internal_energy`, `calc_enthalpy`, `calc_entropy`, `calc_cv`, `calc_cp`, `calc_speed_of_sound`
//! - `calc_joule_momson`, `calc_isothermal_throttling`, `calc_isentropic_temp_pressure`
///
/// Where:
/// - δ = ρ/ρc (reduced density)
/// - τ = Tc/T (inverse reduced temperature)
/// - φ°_τ = ∂φ°/∂τ (second derivative of ideal gas Helmholtz free energy w.r.t. τ)
/// - φʳ_ττ = ∂²φʳ/∂τ² (second derivative of residual Helmholtz free energy w.r.t. τ)
/// - φʳ_δ = ∂φʳ/∂δ (first derivative of residual Helmholtz free energy w.r.t. δ)
/// - φʳ_δδ = ∂²φ/∂δ² (second derivative w.r.t. δ)
/// - φʳ_δτ = ∂²φʳ/∂δτ (mixed second derivative)
/// - φʳ_ττ = ∂²φʳ/∂τ² (second derivative of residual Helmholtz free energy w.r.t. τ)

use crate::iapws95_ideal::*;
use crate::iapws95_residual::*;

// ==========================================================================
// Reference Constants (IAPWS-95 Section 2)
// ==========================================================================

/// Critical temperature: Tc = 647.096 K
pub const IAPWS95_TCRIT: f64 = 647.096;

/// Critical density: rho_c = 322 kg/m³
pub const IAPWS95_RHOCRIT: f64 = 322.0;

/// Critical pressure: pc = 22.064 MPa
pub const IAPWS95_PCRIT: f64 = 22.064;

/// Specific gas constant: R = 0.46151805 kJ/(kg·K)
pub const IAPWS95_R: f64 = 0.46151805;

// ==========================================================================
// Valid Range (IAPWS-95 Section 5)
// ==========================================================================

/// Minimum temperature for practical use (triple point): T_min = 273.16 K
pub const IAPWS95_TMIN: f64 = 273.16;

/// Maximum temperature: Tmax = 1273 K
pub const IAPWS95_TMAX: f64 = 1273.0;

/// Maximum pressure: pmax = 1000 MPa (extended usable range: 100000 MPa)
pub const IAPWS95_PMAX: f64 = 1000.0;

// ==========================================================================
// Helper Functions - Reduced Properties
// ==========================================================================

/// Calculate reduced density delta = rho/rho_c
#[inline]
pub fn reduced_density(rho: f64) -> f64 {
    rho / IAPWS95_RHOCRIT
}

/// Calculate inverse reduced temperature tau = Tc/T
#[inline]
pub fn inv_reduced_temp(T: f64) -> f64 {
    IAPWS95_TCRIT / T
}

// ==========================================================================
// Property Calculations from Helmholtz Free Energy Derivatives
// Based on Table 3 relations of IAPWS-95
// ==========================================================================

/// Compute pressure: p = R·T·ρ·(1 + δ·∂φʳ/∂δ) / 1000 \[MPa\]
#[inline]
pub(crate) fn calc_pressure(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    IAPWS95_R * T * rho * (1.0 + delta * dphi_r_ddelta) / 1000.0
}

/// Compute specific internal energy: u = R·T·τ·(∂φ°/∂τ + ∂φʳ/∂τ) \[kJ/kg\]
#[inline]
pub(crate) fn calc_internal_energy(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    let dphi_dtau = dphi_residual_dtau(delta, tau) + dphi_ideal_dtau(tau);
    IAPWS95_R * T * tau * dphi_dtau
}

/// Compute specific entropy: s = R·[τ·(∂φ°/∂τ + ∂φʳ/∂τ) - φ° - φʳ] \[kJ/(kg·K)\]
#[inline]
pub(crate) fn calc_entropy(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    let phi_o = phi_ideal(delta, tau);
    let phi_r = phi_residual(delta, tau);
    let dphi_dtau = dphi_ideal_dtau(tau)+dphi_residual_dtau(delta, tau);
    IAPWS95_R * (tau * dphi_dtau - phi_o - phi_r)
}

/// Compute specific enthalpy: h = R·T·[τ·(∂φ°/∂τ + ∂φʳ/∂τ) + 1 + δ·∂φʳ/∂δ] \[kJ/kg\]
#[inline]
pub(crate) fn calc_enthalpy(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    let dphi_o_dtau = dphi_ideal_dtau(tau);
    let dphi_r_dtau = dphi_residual_dtau(delta, tau);
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    IAPWS95_R * T * (tau * (dphi_o_dtau + dphi_r_dtau) + 1.0 + delta * dphi_r_ddelta)
}

/// Compute isochoric heat capacity: cv = -R·τ²·(∂²φ°/∂τ² + ∂²φʳ/∂τ²) \[kJ/(kg·K)\]
#[inline]
pub(crate) fn calc_cv(T: f64, rho: f64) -> f64 {
    let tau = inv_reduced_temp(T);
    let delta = reduced_density(rho);
    let d2phi_o_tau = d2phi_ideal_dtau2(tau);
    let d2phi_r_tau2 = d2phi_residual_dtau2(delta, tau);
    IAPWS95_R * (-tau * tau * (d2phi_o_tau + d2phi_r_tau2))
}

/// Compute isobaric heat capacity  kJ/(kg·K):
///  cp/R = -τ²·(∂²φ°/∂τ² + ∂²φʳ/∂τ²) + 
///        (1 + δ*(∂φʳ/∂δ) - δ*τ*(∂²φʳ/∂δ∂τ))² / (1 + 2δ*(∂φʳ/∂δ) + δ²*(∂²φʳ/∂δ²))
#[inline]
pub(crate) fn calc_cp(T: f64, rho: f64) -> f64 {
    let tau = inv_reduced_temp(T);
    let delta = reduced_density(rho);
    let dphi_ddelta = dphi_residual_ddelta(delta, tau);
    let d2phi_ddelta2 = d2phi_residual_ddelta2(delta, tau);
    let d2phi_ddelta_dtau = d2phi_residual_ddelta_dtau(delta, tau);
    let d2phi_o_tau2 = d2phi_ideal_dtau2(tau);
    let d2phi_r_tau2 = d2phi_residual_dtau2(delta, tau);

    // ·τ²·(φ°_ττ + φʳ_ττ)
    let cv_part = -tau * tau * (d2phi_o_tau2 + d2phi_r_tau2);
    //  (1 + δ*φʳ_δ - δ*τ*φʳ_δτ)² / (1 + 2δ*φʳ_δ + δ²*φʳ_δδ)
    let term_temp=1.0 + delta *( dphi_ddelta -  tau * d2phi_ddelta_dtau);
    let numerator =  term_temp* term_temp;
    let denominator = 1.0 + delta *(2.0 *  dphi_ddelta +  delta * d2phi_ddelta2);

    IAPWS95_R * (cv_part + numerator / denominator)
}

/// Compute speed of sound m/s: 
/// w²/RT = 1 + 2δ·∂φʳ/∂δ + δ²·∂²φʳ/∂δ² - N²/(τ²·(∂²φ°/∂τ²+∂²φʳ/∂τ²))
///         N = 1 + δ·∂φʳ/∂δ - δ·τ·∂²φʳ/∂δ∂τ
#[inline]
pub(crate) fn calc_speed_of_sound(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    let d2phi_r_ddelta2 = d2phi_residual_ddelta2(delta, tau);
    let d2phi_r_ddelta_dtau = d2phi_residual_ddelta_dtau(delta, tau);
    let d2phi_o_dtau2 = d2phi_ideal_dtau2(tau);
    let d2phi_r_dtau2 = d2phi_residual_dtau2(delta, tau);
    
    // (1 + δ*φʳ_δ - δ*τ*φʳ_δτ)²
    let term_numerator=1.0 + delta * dphi_r_ddelta - delta * tau * d2phi_r_ddelta_dtau;
    let numerator = term_numerator*term_numerator;
    // (τ²*(φ°_ττ + φʳ_ττ))
    let denominator = tau * tau * (d2phi_o_dtau2 + d2phi_r_dtau2);
    // w² 
     let w_squared = IAPWS95_R * T * (
        1.0 + delta *(2.0 * dphi_r_ddelta +  delta * d2phi_r_ddelta2)
        - numerator / denominator
    );    
    // Convert from kJ/kg to J/kg (multiply by 1000) then take sqrt for m/s
    (w_squared * 1000.0).sqrt()
}

/// Compute Joule-Thomson coefficient: μ = (∂T/∂p)_H [K/MPa]
/// μ·R·ρ = Numerator/Denominator
///     Numerator: -(δ·∂φʳ/∂δ+ δ²·∂²φʳ/∂δ² + δ·∂²φʳ/∂δ∂τ)
///     Denominator: (1+(δ·φʳ/∂δ-δ·τ·∂²φʳ/∂δ∂τ))²
///               -τ²·(∂²φ°/∂τ²+∂²φʳ/∂τ²)·(1 + 2δφʳ_δ + δ²φʳ_δδ)
/// 
#[inline]
pub(crate) fn calc_joule_thomson(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    let d2phi_r_ddelta2 = d2phi_residual_ddelta2(delta, tau);
    let d2phi_r_ddelta_dtau = d2phi_residual_ddelta_dtau(delta, tau);
    let d2phi_o_dtau2 = d2phi_ideal_dtau2(tau);
    let d2phi_r_dtau2 = d2phi_residual_dtau2(delta, tau);
    // Numerator: 
    //  -(δ·∂φʳ/∂δ+ δ²·∂²φʳ/∂δ² + δ·∂²φʳ/∂δ∂τ)
    let numerator = -delta * (dphi_r_ddelta + delta*d2phi_r_ddelta2 + d2phi_r_ddelta_dtau);
    // Denominator: 
    // (1+(δ·φʳ/∂δ-δ·τ·∂²φʳ/∂δ∂τ))²
    // -τ²·(∂²φ°/∂τ²+∂²φʳ/∂τ²)·(1 + 2δφʳ_δ + δ²φʳ_δδ)
    //
    let term_1 = 1.0 + delta * (dphi_r_ddelta -  tau * d2phi_r_ddelta_dtau);
    let term_1_2 = term_1*term_1;
    // τ²·(∂²φ°/∂τ²+∂²φʳ/∂τ²)·(1 + 2δφʳ_δ + δ²φʳ_δδ)
    let term2 = tau*tau*(d2phi_o_dtau2+d2phi_r_dtau2)
                     *(1.0 +  delta *(2.0 * dphi_r_ddelta +  delta * d2phi_r_ddelta2));
    let denominator = IAPWS95_R *rho *(term_1_2- term2);
    
    1000.0*(numerator / denominator)

}

/// Compute Isothermal throttling coefficient: (∂τ/∂p)_T kJ/(kg·MPa)
/// 
/// Based on IAPWS-95 Table 3 relations:
/// (∂τ/∂p)_T = 1-(1 + δφʳ_δ - δτφʳ_δτ) / (1 + 2δφʳ_δ + δ²φʳ_δδ)
/// 
/// Where:
/// - δ = ρ/ρc (reduced density)
/// - τ = Tc/T (inverse reduced temperature)
/// - φʳ_δ = ∂φʳ/∂δ (first derivative of residual Helmholtz free energy w.r.t. δ)
/// - φʳ_δδ = ∂²φʳ/∂δ² (second derivative w.r.t. δ)
/// - φʳ_δτ = ∂²φʳ/∂δ∂τ (mixed second derivative)
#[inline]
pub(crate) fn calc_isothermal_throttling(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    let d2phi_r_ddelta2 = d2phi_residual_ddelta2(delta, tau);
    let d2phi_r_ddelta_dtau = d2phi_residual_ddelta_dtau(delta, tau);
    
    // Numerator: 1 + δφʳ_δ - δτφʳ_δτ
    let numerator = 1.0 + delta * (dphi_r_ddelta -  tau * d2phi_r_ddelta_dtau);
    
    // Denominator: 1 + 2δφʳ_δ + δ²φʳ_δδ
    let denominator = 1.0 + 2.0 * delta *( dphi_r_ddelta +  delta * d2phi_r_ddelta2);
    
    // (∂τ/p)_T = 1 - (1 + δφʳ_δ - δτφʳ_δτ) / (1 + 2δφ_δ + δ²φʳ_δδ)
    1.0 - (numerator / denominator)
}

/// Compute Isentropic temperature-pressure coefficient: β_s = (∂T/∂p)_s 1/K
/// 
/// Based on IAPWS-95 Table 3 relations:
/// β_s * ρ * R = (1 + δφ_δ - δτφʳ_δτ) / [(1 + δφʳ_δ - δτφʳ_δτ)² - τ²(φ°_ττ + φʳ_ττ)(1 + 2δφʳ_δ + δ²φʳ_δδ)]
/// 
#[inline]
pub(crate) fn calc_isentropic_temp_pressure(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    let d2phi_r_ddelta2 = d2phi_residual_ddelta2(delta, tau);
    let d2phi_r_ddelta_dtau = d2phi_residual_ddelta_dtau(delta, tau);
    let d2phi_o_dtau2 = d2phi_ideal_dtau2(tau);
    let d2phi_r_dtau2 = d2phi_residual_dtau2(delta, tau);
    
    // Numerator: 1 + δφʳ_δ - δτφʳ_δτ
    let numerator = 1.0 + delta * dphi_r_ddelta - delta * tau * d2phi_r_ddelta_dtau;
    
    // Denominator: (1 + δφʳ_δ - δτφʳ_δτ)² - τ²(φ°_ττ + φʳ_ττ)(1 + 2δφʳ_δ + δ²φʳ_δδ)
    let term1 = 1.0 + delta * dphi_r_ddelta - delta * tau * d2phi_r_ddelta_dtau;
    let term2 = 1.0 + 2.0 * delta * dphi_r_ddelta + delta * delta * d2phi_r_ddelta2;
    let denominator = (term1 * term1) - tau * tau * (d2phi_o_dtau2 + d2phi_r_dtau2) * term2;
    
    // β_s = numerator / (ρ * R * denominator)
    numerator / (rho * IAPWS95_R * denominator)
}

// ==========================================================================
// Main API Functions
// ==========================================================================

/// Check if a state is within the valid range.
#[inline]
pub fn iapws95_in_range(T: f64, _p: Option<f64>) -> bool {
    let _ = _p;
    T >= IAPWS95_TMIN && T <= IAPWS95_TMAX
}

// ==========================================================================
// Convenience Functions for (t_c,rho) → property calculations with °C input
// ==========================================================================

/// Calculate pressure at given temperature (°C) and density \[MPa\]
#[inline]
pub fn tr2p(t_c: f64, rho: f64) -> f64 {
    let t_k = t_c + 273.15;
    calc_pressure(t_k, rho)
}

/// Calculate internal energy at given temperature (°C) and density \[kJ/kg\]
#[inline]
pub fn tr2u(t_c: f64, rho: f64) -> f64 {
    let t_k = t_c + 273.15;
    calc_internal_energy(t_k, rho)
}

/// Calculate enthalpy at given temperature (°C) and density \[kJ/kg\]
#[inline]
pub fn tr2h(t_c: f64, rho: f64) -> f64 {
    let t_k = t_c + 273.15;
    calc_enthalpy(t_k, rho)
}

/// Calculate entropy at given temperature (°C) and density \[kJ/(kg·K)\]
#[inline]
pub fn tr2s(t_c: f64, rho: f64) -> f64 {
    let t_k = t_c + 273.15;
    calc_entropy(t_k, rho)
}

/// Calculate constant-volume specific heat at given temperature (°C) and density \[kJ/(kg·K)\]
#[inline]
pub fn tr2cv(t_c: f64, rho: f64) -> f64 {
    let t_k = t_c + 273.15;
    calc_cv(t_k, rho)
}

/// Calculate constant-pressure specific heat at given temperature (°C) and density \[kJ/(kg·K)\]
#[inline]
pub fn tr2cp(t_c: f64, rho: f64) -> f64 {
    let t_k = t_c + 273.15;
    calc_cp(t_k, rho)
}

/// Calculate speed of sound at given temperature (°C) and density \[m/s\]
#[inline]
pub fn tr2w(t_c: f64, rho: f64) -> f64 {
    let t_k = t_c + 273.15;
    calc_speed_of_sound(t_k, rho)
}

/// Calculate Joule-Thomson coefficient at given temperature (°C) and density \[K/MPa\]
#[inline]
pub fn tr2jt(t_c: f64, rho: f64) -> f64 {
    let t_k = t_c + 273.15;
    calc_joule_thomson(t_k, rho)
}

/// Calculate Isothermal throttling coefficient at given temperature (°C) and density \[kJ/(kg·MPa)\]
#[inline]
pub fn tr2itt(t_c: f64, rho: f64) -> f64 {
    let t_k = t_c + 273.15;
    calc_isothermal_throttling(t_k, rho)
}

/// Calculate Isentropic temperature-pressure coefficient at given temperature (°C) and density \1/K]
#[inline]
pub fn tr2beta_s(t_c: f64, rho: f64) -> f64 {
    let t_k = t_c + 273.15;
    calc_isentropic_temp_pressure(t_k, rho)
}