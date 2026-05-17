//! IAPWS-95 Header - Reference Constants, Ranges, and Data Structures
//!
//! Translated from iapws95.h
use crate::iapws95_ideal::*;
use crate::iapws95_residual::*;

// ==========================================================================
// Reference Constants (IAPWS-95 Section 2)
// ==========================================================================

/// Critical temperature: Tc = 647.096 K
pub const IAPWS95_TCRIT: f64 = 647.096;

/// Critical density: rho_c = 322 kg/m³
pub const IAPWS95_RHOCRIT: f64 = 322.0;

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

/// Compute pressure: p = RT*delta*(1 + delta*ddelta) [MPa]
pub fn calc_pressure(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    IAPWS95_R * T * rho * (1.0 + delta * dphi_r_ddelta) / 1000.0
}

/// Compute specific internal energy: u = RT*tau*(dphi_o/dtau + dphi_r/dtau) [kJ/kg]
pub fn calc_internal_energy(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    let dphi_dtau = dphi_residual_dtau(delta, tau) + dphi_ideal_dtau(tau);
    IAPWS95_R * T * tau * dphi_dtau
}

/// Compute specific entropy: s = R*(tau*dphi/dtau - phi_o - phi_r) [kJ/(kg*K)]
pub fn calc_entropy(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    let phi_o = phi_ideal(delta, tau);
    let phi_r = phi_residual(delta, tau);
    let phi_o_t = dphi_ideal_dtau(tau);
    let phi_r_t = dphi_residual_dtau(delta, tau);
    let dphi_dtau = phi_o_t + phi_r_t;
    IAPWS95_R * (tau * dphi_dtau - phi_o - phi_r)
}

/// Compute specific enthalpy: h = RT*[tau*(dphi_o/dtau + dphi_r/dtau) + 1 + delta*(1 + delta*dphi_r/ddelta)] [kJ/kg]
pub fn calc_enthalpy(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    let dphi_o_dtau = dphi_ideal_dtau(tau);
    let dphi_r_dtau = dphi_residual_dtau(delta, tau);
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    IAPWS95_R * T * (tau * (dphi_o_dtau + dphi_r_dtau) + 1.0 + delta * dphi_r_ddelta)
}

/// Compute isochoric heat capacity: cv = R*(-tau^2*(d2phi_o_tau2+d2phi_r_tau2)) [kJ/(kg*K)]
pub fn calc_cv(T: f64, rho: f64) -> f64 {
    let tau = inv_reduced_temp(T);
    let delta = reduced_density(rho);
    let phi_o_tt = d2phi_ideal_dtau2(tau);
    let phi_r_tt = d2phi_residual_dtau2(delta, tau);
    IAPWS95_R * (-tau * tau * (phi_o_tt + phi_r_tt))
}

/// Compute isobaric heat capacity: cp = cv + R*(1 + δ*(∂φʳ/∂δ) - δ*τ*(∂²φʳ/∂δ∂τ))² / (1 + 2δ*(∂φʳ/∂δ) + δ²*(∂²φʳ/∂δ²)) [kJ/(kg*K)]
pub fn calc_cp(T: f64, rho: f64) -> f64 {
    let tau = inv_reduced_temp(T);
    let delta = reduced_density(rho);
    let dphi_ddelta = dphi_residual_ddelta(delta, tau);
    let d2phi_ddelta2 = d2phi_residual_ddelta2(delta, tau);
    let d2phi_ddelta_dtau = d2phi_residual_ddelta_dtau(delta, tau);
    
    let cv_val = calc_cv(T, rho);
    
    // cp = cv + R * (1 + δ*φʳ_δ - δ*τ*φʳ_δτ)² / (1 + 2δ*φʳ_δ + δ²*φʳ_δδ)
    let numerator = (1.0 + delta * dphi_ddelta - delta * tau * d2phi_ddelta_dtau).powi(2);
    let denominator = 1.0 + 2.0 * delta * dphi_ddelta + delta * delta * d2phi_ddelta2;
    
    cv_val + IAPWS95_R * numerator / denominator
}

/// Compute speed of sound: w [m/s]
pub fn calc_speed_of_sound(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    
    let dphi_ddelta = dphi_residual_ddelta(delta, tau);
    let d2phi_ddelta2 = d2phi_residual_ddelta2(delta, tau);
    let d2phi_ddelta_dtau = d2phi_residual_ddelta_dtau(delta, tau);
    let d2phi_dtau2_ideal = d2phi_ideal_dtau2(tau);
    let d2phi_dtau2_residual = d2phi_residual_dtau2(delta, tau);
    
    // w² = R*T * [1 + 2δ*φʳ_δ + δ²*φʳ_δδ - (1 + δ*φʳ_δ - δ*τ*φʳ_δτ)² / (τ²*(φ°_ττ + φʳ_ττ))]
    let numerator = (1.0 + delta * dphi_ddelta - delta * tau * d2phi_ddelta_dtau).powi(2);
    let denominator = tau * tau * (d2phi_dtau2_ideal + d2phi_dtau2_residual);
    
    let w_squared = IAPWS95_R * T * (
        1.0 + 2.0 * delta * dphi_ddelta + delta * delta * d2phi_ddelta2 
        - numerator / denominator
    );
    
    // Convert from kJ/kg to J/kg (multiply by 1000) then take sqrt for m/s
    (w_squared * 1000.0).sqrt()
}

// ==========================================================================
// Main API Functions
// ==========================================================================

/// Check if a state is within the valid range.
pub fn iapws95_in_range(T: f64, _p: Option<f64>) -> bool {
    let _ = _p;
    T >= IAPWS95_TMIN && T <= IAPWS95_TMAX
}