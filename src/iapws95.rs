//! IAPWS-95 Header - Reference Constants, Ranges, and Data Structures
//!
//! Translated from iapws95.h
use  crate::iapws95_ideal::*;
use  crate::iapws95_residual::*;

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
    let dphi_r_ddelta = dphi_residual_ddelta(delta, inv_reduced_temp(T));
    IAPWS95_R * T * delta * (1.0 + delta * dphi_r_ddelta)
}

/// Compute specific internal energy: u = RT*tau*(phi_o + phi_r + tau*dphi/dtau) [kJ/kg]
pub fn calc_internal_energy(T: f64, rho: f64) -> f64 {
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);

    let phi_o = phi_ideal(delta, tau);
    let phi_r =phi_residual(delta, tau);
    let dphi_dtau =
        dphi_residual_dtau(delta, tau) + crate::iapws95_ideal::dphi_ideal_dtau(tau);
    IAPWS95_R * T * tau * (phi_o + phi_r + tau * dphi_dtau)
}

/// Compute specific entropy: s = R*(phi_o + phi_r - tau*dphi/dtau) [kJ/(kg*K)]
pub fn calc_entropy(T: f64, rho:f64) -> f64{
    let delta = reduced_density(rho);
    let tau = inv_reduced_temp(T);
    let phi_o = phi_ideal(delta, tau);
    let phi_r = phi_residual(delta, tau);
    let phi_o_t = dphi_ideal_dtau(tau);
    let phi_r_t = dphi_residual_dtau(delta, tau);
    let dphi_dtau = phi_o_t + phi_r_t;
    IAPWS95_R * (phi_o + phi_r - dphi_dtau)
}

/// Compute specific enthalpy: h = u + p/rho [kJ/kg]
pub fn calc_enthalpy(_T: f64, rho: f64, p: f64, u: f64) -> f64 {
    u + p / rho * 1000.0 // p in MPa, rho in kg/m3: p/rho*1000 converts to kJ/kg
}

/// Compute isochoric heat capacity: cv = R*(-tau^2*d2phi/dtau2) [kJ/(kg*K)]
pub fn calc_cv(T: f64, d2phi_dtau2: f64) -> f64 {
    IAPWS95_R * (-T * T * d2phi_dtau2)
}

/// Compute isobaric heat capacity: cp = cv + R*(dp/dT)_rho^2 / (dp/drho)_T [kJ/(kg*K)]
pub fn calc_cp(T: f64, rho: f64, dphi_ddelta: f64, d2phi_ddelta2: f64) -> f64 {
    let delta = reduced_density(rho);

    // Partial derivatives needed for cp calculation
    let dpdT = IAPWS95_R * (1.0 + dphi_ddelta); // (dp/dT)_rho in MPa/K
    let dpdrho = IAPWS95_R * T * delta / rho * (3.0 + d2phi_ddelta2); // (dp/drho)_T in MPa/(kg/m3)

    let cv_val = calc_cv(T, -1.0); // Placeholder - would need actual d2phi/dtau2

    cv_val + IAPWS95_R * dpdT * dpdT / dpdrho
}

/// Compute speed of sound: w [m/s]
pub fn calc_speed_of_sound(_rho: f64, T: f64, cp: f64, cv: f64) -> f64 {
    // Simplified formula using heat capacity ratio and basic compressibility
    (cp / cv * IAPWS95_R * T * 1000.0).sqrt()
}

// ==========================================================================
// Main API Functions
// ==========================================================================

/// Check if a state is within the valid range.
pub fn iapws95_in_range(T: f64, _p: Option<f64>) -> bool {
    let _ = _p;
    T >= IAPWS95_TMIN && T <= IAPWS95_TMAX
}