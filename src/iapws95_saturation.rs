//! IAPWS-95 Saturation Properties Calculation
//!
//! Computes saturation properties using a hybrid approach:
//! 1. Use IAPWS SR1-86 (1992) explicit equations for initial guesses
//! 2. Refine using Newton's method to solve IAPWS-95 phase equilibrium conditions
//!
//! Phase equilibrium conditions:
//! - Equal pressure: p(δ', τ) = p(δ'', τ)
//! - Equal chemical potential: μ(δ', τ) = μ(δ'', τ)
//!
//! This approach combines the robustness of SR1-86 initial estimates with
//! the accuracy of IAPWS-95.

use crate::iapws95::*;
use crate::iapws95_residual::*;

/// SR1-86 reference constants for initial guesses
const SR1_TC: f64 = 647.096;
const SR1_RHOC: f64 = 322.0;

/// Saturation properties at a given temperature
pub struct SaturationProperties {
    /// Saturation vapor pressure [MPa]
    pub p_sat: f64,
    /// Saturated liquid density [kg/m³]
    pub rho_l: f64,
    /// Saturated vapor density [kg/m³]
    pub rho_v: f64,
    /// Saturated liquid specific enthalpy [kJ/kg]
    pub h_l: f64,
    /// Saturated vapor specific enthalpy [kJ/kg]
    pub h_v: f64,
    /// Saturated liquid specific entropy [kJ/(kg·K)]
    pub s_l: f64,
    /// Saturated vapor specific entropy [kJ/(kg·K)]
    pub s_v: f64,
}

/// SR1-86: Compute saturated liquid density for initial guess
fn sr1_liquid_density(t: f64) -> f64 {
    let theta = t / SR1_TC;
    let tau = 1.0 - theta;
    
    let b1 = 1.99274064;
    let b2 = 1.09965342;
    let b3 = -0.510839303;
    let b4 = -1.75493479;
    let b5 = -45.5170352;
    let b6 = -6.74694450e5;
    
    let rho_ratio = 1.0
        + b1 * tau.powf(1.0 / 3.0)
        + b2 * tau.powf(2.0 / 3.0)
        + b3 * tau.powf(5.0 / 3.0)
        + b4 * tau.powf(16.0 / 3.0)
        + b5 * tau.powf(43.0 / 3.0)
        + b6 * tau.powf(110.0 / 3.0);
    
    SR1_RHOC * rho_ratio
}

/// SR1-86: Compute saturated vapor density for initial guess
fn sr1_vapor_density(t: f64) -> f64 {
    let theta = t / SR1_TC;
    let tau = 1.0 - theta;
    
    let c1 = -2.03150240;
    let c2 = -2.68302940;
    let c3 = -5.38626492;
    let c4 = -17.2991605;
    let c5 = -44.7586581;
    let c6 = -63.9201063;
    
    let ln_rho_ratio = c1 * tau.powf(2.0 / 6.0)
        + c2 * tau.powf(4.0 / 6.0)
        + c3 * tau.powf(8.0 / 6.0)
        + c4 * tau.powf(18.0 / 6.0)
        + c5 * tau.powf(37.0 / 6.0)
        + c6 * tau.powf(71.0 / 6.0);
    
    SR1_RHOC * ln_rho_ratio.exp()
}

/// Compute dimensionless pressure term: J = δ·(1 + δ·∂φʳ/∂δ)
#[inline]
fn pressure_term(delta: f64, tau: f64) -> f64 {
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    delta * (1.0 + delta * dphi_r_ddelta)
}

/// Compute dimensionless chemical potential term: K = δ·∂φʳ/∂δ + φʳ + ln(δ)
#[inline]
fn chemical_potential_term(delta: f64, tau: f64) -> f64 {
    let phi_r = phi_residual(delta, tau);
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    delta * dphi_r_ddelta + phi_r + delta.ln()
}

/// Compute derivatives (dJ/dδ, dK/dδ)
fn compute_derivatives(delta: f64, tau: f64) -> (f64, f64) {
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    let d2phi_r_ddelta2 = d2phi_residual_ddelta2(delta, tau);
    
    let dj_ddelta = 1.0 + 2.0 * delta * dphi_r_ddelta + delta * delta * d2phi_r_ddelta2;
    let dk_ddelta = 2.0 * dphi_r_ddelta + delta * d2phi_r_ddelta2 + 1.0 / delta;
    
    (dj_ddelta, dk_ddelta)
}

/// Solve phase equilibrium using Newton's method with SR1-86 initial guesses
/// 
/// Solves:
///   F1(δL, δV) = K(δV, τ) - K(δL, τ) = 0  (equal chemical potential)
///   F2(δL, δV) = J(δV, τ) - J(δL, τ) = 0  (equal pressure)
fn solve_phase_equilibrium(t: f64, tau: f64) -> Option<(f64, f64)> {
    // Use SR1-86 equations for initial guesses
    let rho_l_init = sr1_liquid_density(t);
    let rho_v_init = sr1_vapor_density(t);
    
    let mut delta_l = rho_l_init / IAPWS95_RHOCRIT;
    let mut delta_v = rho_v_init / IAPWS95_RHOCRIT;
    
    // Clamp to reasonable ranges
    delta_l = delta_l.clamp(1.0, 4.5);
    delta_v = delta_v.clamp(1e-7, 0.8);
    
    // Newton iteration with damping
    for iter in 0..100 {
        let jl = pressure_term(delta_l, tau);
        let jv = pressure_term(delta_v, tau);
        let kl = chemical_potential_term(delta_l, tau);
        let kv = chemical_potential_term(delta_v, tau);
        
        let f1 = kv - kl;
        let f2 = jv - jl;
        
        if f1.abs() < 1e-12 && f2.abs() < 1e-12 {
            return Some((delta_l, delta_v));
        }
        
        let (dj_dl, dk_dl) = compute_derivatives(delta_l, tau);
        let (dj_dv, dk_dv) = compute_derivatives(delta_v, tau);
        
        let j11 = -dk_dl;
        let j12 = dk_dv;
        let j21 = -dj_dl;
        let j22 = dj_dv;
        
        let det = j11 * j22 - j12 * j21;
        
        if det.abs() < 1e-15 {
            break;
        }
        
        let delta_delta_l = -(j22 * f1 - j12 * f2) / det;
        let delta_delta_v = -(j11 * f2 - j21 * f1) / det;
        
        let mut damping = 1.0;
        let mut new_delta_l = delta_l + damping * delta_delta_l;
        let mut new_delta_v = delta_v + damping * delta_delta_v;
        
        for _ in 0..20 {
            if new_delta_l > 0.8 && new_delta_l < 5.0 
                && new_delta_v > 1e-8 && new_delta_v < 2.0
                && new_delta_l > new_delta_v * 1.1
            {
                break;
            }
            damping *= 0.5;
            new_delta_l = delta_l + damping * delta_delta_l;
            new_delta_v = delta_v + damping * delta_delta_v;
        }
        
        if new_delta_l <= 0.8 || new_delta_l >= 5.0 
            || new_delta_v <= 1e-8 || new_delta_v >= 2.0
            || new_delta_l <= new_delta_v * 1.1
        {
            break;
        }
        
        delta_l = new_delta_l;
        delta_v = new_delta_v;
        
        if iter > 5 {
            let new_jl = pressure_term(delta_l, tau);
            let new_jv = pressure_term(delta_v, tau);
            let new_kl = chemical_potential_term(delta_l, tau);
            let new_kv = chemical_potential_term(delta_v, tau);
            
            if (new_kv - new_kl).abs() < 1e-10 && (new_jv - new_jl).abs() < 1e-10 {
                return Some((delta_l, delta_v));
            }
        }
    }
    
    Some((delta_l, delta_v))
}

/// Compute saturation properties at given temperature T: [K] 
pub fn calc_saturation_properties(T: f64) -> Option<SaturationProperties> {
    if T < 273.16 || T > IAPWS95_TCRIT {
        return None;
    }

    let tau = inv_reduced_temp(T);
    let (delta_l, delta_v) = solve_phase_equilibrium(T, tau)?;

    let rho_l = delta_l * IAPWS95_RHOCRIT;
    let rho_v = delta_v * IAPWS95_RHOCRIT;

    let p_sat = calc_pressure(T, rho_l);
    let h_l = calc_enthalpy(T, rho_l);
    let h_v = calc_enthalpy(T, rho_v);
    let s_l = calc_entropy(T, rho_l);
    let s_v = calc_entropy(T, rho_v);

    Some(SaturationProperties {
        p_sat,
        rho_l,
        rho_v,
        h_l,
        h_v,
        s_l,
        s_v,
    })
}

/// Compute saturation properties at given temperature t_c,°C 
pub fn sat_t(t_c: f64) -> Option<SaturationProperties> {
    if t_c < 0.01 || t_c > IAPWS95_TCRIT - 273.15 {
        return None;
    }

    let t_k = t_c + 273.15;
    calc_saturation_properties(t_k)
}