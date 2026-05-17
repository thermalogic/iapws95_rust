//! IAPWS-95 Inverse Problem Solver - (p,T) → ρ
//! 
//! Provides numerical solution for density when pressure and temperature are known.
//! This is essential for applications where (p,T) conditions are given instead of (T,ρ).
//! 
//! # Algorithm
//! 
//! Uses Newton's method with multiple initial guesses based on ideal gas law:
//! 1. **Initial guess**: ρ₀ = p/(R·T) from ideal gas equation
//! 2. **Guess selection**: Choose appropriate starting point based on density regime (vapor/liquid region)
//! 3. **Newton iteration**: Solve f(ρ) = p_calc(T,ρ) - p = 0 with adaptive damping

use crate::iapws95::{calc_pressure, IAPWS95_R, IAPWS95_RHOCRIT};

// ==========================================================================
// Public API Functions
// ==========================================================================

/// Solve for density rho at given pressure p [MPa] and temperature T [K]
/// 
/// Uses Newton's method with multiple initial guesses based on ideal gas law.
/// Returns None if convergence fails.
/// 
/// # Arguments
/// * `p` - Pressure in MPa
/// * `T` - Temperature in Kelvin
/// 
/// # Returns
/// * `Some(rho)` - Density in kg/m³ if successful
/// * `None` - If Newton's method fails to converge within 200 iterations
/// 
/// # Examples
/// ```
/// use iapws95::iapws95_pT::solve_density;
/// 
/// // Given p=16.10 MPa, T=808.25 K (535.10°C)
/// if let Some(rho) = solve_density(16.10, 808.25) {
///     println!("Density: {:.4} kg/m³", rho);
/// }
/// ```
/// 
pub fn solve_density(p: f64, T: f64) -> Option<f64> {
    if p <= 0.0 || T <= 0.0 {
        return None;
    }

    // Initial guess from ideal gas law: rho = p*1000/(R*T)
    let rho_ideal = p * 1000.0 / (IAPWS95_R * T);
    
    // Choose initial guesses based on density regime
    let guesses = if rho_ideal < IAPWS95_RHOCRIT * 0.5 {
        // Low density (vapor region)
        vec![rho_ideal.max(1e-6), rho_ideal * 0.5, rho_ideal * 2.0, IAPWS95_RHOCRIT]
    } else {
        // High density (liquid region)
        vec![IAPWS95_RHOCRIT, rho_ideal, rho_ideal * 0.5, rho_ideal * 2.0]
    };

    for &rho_init in &guesses {
        if let Some(rho) = solve_density_newton(p, T, rho_init) {
            return Some(rho);
        }
    }

    None
}

// ==========================================================================
// Private Helper Functions
// ==========================================================================

/// Newton's method iteration for density solver
fn solve_density_newton(p: f64, T: f64, rho_init: f64) -> Option<f64> {
    let mut rho = rho_init.max(1e-8);

    for _ in 0..200 {
        let p_calc = calc_pressure(T, rho);
        let f = p_calc - p;
        
        if f.abs() < 1e-6 {
            return Some(rho);
        }

        let drho = (rho * 1e-8).max(1e-10);
        let p_plus = calc_pressure(T, rho + drho);
        let df_drho = (p_plus - p_calc) / drho;

        if df_drho.abs() < 1e-20 {
            break;
        }

        let delta_rho = -f / df_drho;
        let damping = if delta_rho.abs() > rho * 0.5 { 0.1 } else { 0.5 };
        rho = (rho + damping * delta_rho).max(1e-8);
    }

    let p_calc = calc_pressure(T, rho);
    if (p_calc - p).abs() < 1e-3 {
        Some(rho)
    } else {
        None
    }
}
