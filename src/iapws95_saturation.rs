//! IAPWS-95 Saturation Properties Calculation
//!
//! Computes saturation properties using phase-equilibrium condition:
//! - Equal pressure: p(δ', τ) = p(δ'', τ) = p_σ
//! - Equal Gibbs energy: g(δ', τ) = g(δ'', τ)
//!
//! Based on IAPWS-95 Table 8 and Section 4.

use crate::iapws95::*;
use crate::iapws95_ideal::*;
use crate::iapws95_residual::*;

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

/// Compute dimensionless Helmholtz free energy φ = φ° + φʳ
#[inline]
fn phi_total(delta: f64, tau: f64) -> f64 {
    phi_ideal(delta, tau) + phi_residual(delta, tau)
}

/// Compute dimensionless Gibbs free energy: g/(RT) = φ + δ·(∂φ/∂δ)
#[inline]
fn gibbs_reduced(delta: f64, tau: f64) -> f64 {
    let phi = phi_total(delta, tau);
    let dphi_ddelta = 1.0 / delta + dphi_residual_ddelta(delta, tau);
    phi + delta * dphi_ddelta
}

/// Compute pressure: p = ρ·R·T·(1 + δ·∂φʳ/∂δ) [MPa]
#[inline]
fn calc_pressure_from_delta(delta: f64, tau: f64) -> f64 {
    let t = IAPWS95_TCRIT / tau;
    let dphi_r_ddelta = dphi_residual_ddelta(delta, tau);
    let rho = delta * IAPWS95_RHOCRIT;
    IAPWS95_R * t * rho * (1.0 + delta * dphi_r_ddelta) / 1000.0
}

/// Estimate saturation pressure using Wagner equation
fn estimate_saturation_pressure(tau: f64) -> f64 {
    let t_reduced = 1.0 / tau;
    let theta = 1.0 - t_reduced;
    
    let a1 = -7.85951783;
    let a2 = 1.84408259;
    let a3 = -11.7866497;
    let a4 = 22.6807411;
    let a5 = -15.9618719;
    let a6 = 1.80122502;
    
    let ln_pr = (1.0 / t_reduced) * (
        a1 * theta
        + a2 * theta.powf(1.5)
        + a3 * theta.powf(3.0)
        + a4 * theta.powf(3.5)
        + a5 * theta.powf(4.0)
        + a6 * theta.powf(7.5)
    );
    
    IAPWS95_PCRIT * ln_pr.exp()
}

/// Find density root for given pressure using bisection method
/// Returns delta such that p(delta, tau) = p_target
fn find_density_bisection(p_target: f64, tau: f64, delta_low: f64, delta_high: f64) -> Option<f64> {
    let mut lo = delta_low;
    let mut hi = delta_high;
    
    let p_lo = calc_pressure_from_delta(lo, tau);
    let p_hi = calc_pressure_from_delta(hi, tau);
    
    // Check if target is within range
    if (p_target - p_lo) * (p_target - p_hi) > 0.0 {
        return None;
    }
    
    for _ in 0..300 {
        let mid = (lo + hi) / 2.0;
        let p_mid = calc_pressure_from_delta(mid, tau);
        
        if (p_mid - p_target).abs() < 1e-12 {
            return Some(mid);
        }
        
        if (hi - lo) < 1e-14 {
            return Some(mid);
        }
        
        if p_mid > p_target {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    
    Some((lo + hi) / 2.0)
}

/// Find vapor density root by scanning for monotonic region
fn find_vapor_density(p_target: f64, tau: f64) -> Option<f64> {
    // Scan to find the peak pressure point in vapor region
    let mut p_max = 0.0;
    let mut delta_at_max = 0.01;
    
    for i in 1..=1000 {
        let delta = i as f64 / 1000.0;
        let p = calc_pressure_from_delta(delta, tau);
        if p > p_max {
            p_max = p;
            delta_at_max = delta;
        }
    }
    
    // Search in the monotonic increasing region [1e-6, delta_at_max]
    // But first check if p_target is achievable
    let p_at_max = p_max;
    if p_target > p_at_max {
        return None;
    }
    
    // Find the exact range where pressure crosses p_target
    let mut delta_start = 1e-6;
    let mut delta_end = delta_at_max;
    
    // Binary search for the right range
    for _ in 0..50 {
        let mid = (delta_start + delta_end) / 2.0;
        let p_mid = calc_pressure_from_delta(mid, tau);
        
        if p_mid > p_target {
            delta_end = mid;
        } else {
            delta_start = mid;
        }
    }
    
    // Now search in [delta_start * 0.5, delta_end]
    let lo = delta_start * 0.5;
    let hi = delta_end.min(delta_at_max * 0.99);
    
    find_density_bisection(p_target, tau, lo, hi)
}

/// Find liquid density root by scanning for appropriate region
fn find_liquid_density(p_target: f64, tau: f64) -> Option<f64> {
    // Scan from delta=1.0 to 4.0 to find where pressure crosses p_target
    let mut delta_cross: Option<f64> = None;
    let mut p_prev = calc_pressure_from_delta(1.0, tau);
    
    for i in 100..=400 {
        let delta = i as f64 / 100.0;
        let p = calc_pressure_from_delta(delta, tau);
        
        // Check if we crossed p_target
        if (p_target - p_prev) * (p_target - p) < 0.0 {
            delta_cross = Some(delta - 0.01);
            break;
        }
        p_prev = p;
    }
    
    if let Some(delta_start) = delta_cross {
        find_density_bisection(p_target, tau, delta_start, 4.0)
    } else {
        // Try full range
        find_density_bisection(p_target, tau, 1.0, 4.0)
    }
}

/// Solve phase equilibrium using pressure bisection
/// Returns (delta_liquid, delta_vapor) for given tau
fn solve_phase_equilibrium(tau: f64) -> Option<(f64, f64)> {
    // Initial pressure estimate from Wagner equation
    let p_est = estimate_saturation_pressure(tau);
    
    if p_est <= 0.0 || p_est > IAPWS95_PCRIT {
        return None;
    }
    
    let mut best_delta_l: Option<f64> = None;
    let mut best_delta_v: Option<f64> = None;
    let mut best_dg_abs = f64::MAX;
    
    // Iterative refinement using pressure bisection
    let mut p_low = p_est * 0.5;
    let mut p_high = p_est * 2.0;
    
    for _ in 0..2000 {
        let p_mid = (p_low + p_high) / 2.0;
        
        // Find density roots using dynamic range detection
        let delta_l = find_liquid_density(p_mid, tau);
        let delta_v = find_vapor_density(p_mid, tau);
        
        if delta_l.is_none() || delta_v.is_none() {
            // Adjust pressure bounds
            if delta_l.is_none() {
                p_high = p_mid;
            } else {
                p_low = p_mid;
            }
            continue;
        }
        
        let delta_l = delta_l.unwrap();
        let delta_v = delta_v.unwrap();
        
        if delta_l <= delta_v {
            p_high = p_mid;
            continue;
        }
        
        // Check phase equilibrium: g_l = g_v
        let g_l = gibbs_reduced(delta_l, tau);
        let g_v = gibbs_reduced(delta_v, tau);
        let dg = g_l - g_v;
        let dg_abs = dg.abs();
        
        if dg_abs < best_dg_abs {
            best_dg_abs = dg_abs;
            best_delta_l = Some(delta_l);
            best_delta_v = Some(delta_v);
        }
        
        if dg_abs < 1e-10 {
            return Some((delta_l, delta_v));
        }
        
        // Adjust pressure based on Gibbs energy difference
        if dg > 0.0 {
            p_low = p_mid;
        } else {
            p_high = p_mid;
        }
    }
    
    if let (Some(dl), Some(dv)) = (best_delta_l, best_delta_v) {
        if dl > dv {
            Some((dl, dv))
        } else {
            None
        }
    } else {
        None
    }
}

/// Compute saturation properties at given temperature
/// 
/// # Arguments
/// * `T` - Temperature [K], must be in range [273.16, 647.096]
/// 
/// # Returns
/// * `Some(SaturationProperties)` if calculation succeeds
/// * `None` if T is out of valid saturation range
pub fn calc_saturation_properties(t: f64) -> Option<SaturationProperties> {
    if t < IAPWS95_TMIN || t > IAPWS95_TCRIT {
        return None;
    }

    let tau = inv_reduced_temp(t);
    let (delta_l, delta_v) = solve_phase_equilibrium(tau)?;

    let rho_l = delta_l * IAPWS95_RHOCRIT;
    let rho_v = delta_v * IAPWS95_RHOCRIT;

    let p_sat = calc_pressure(t, rho_l);
    let h_l = calc_enthalpy(t, rho_l);
    let h_v = calc_enthalpy(t, rho_v);
    let s_l = calc_entropy(t, rho_l);
    let s_v = calc_entropy(t, rho_v);

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
