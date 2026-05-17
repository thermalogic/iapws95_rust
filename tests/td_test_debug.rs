#![allow(warnings)]
///  cargo test --test td_test_debug -- --nocapture
use iapws95::iapws95::*;
use iapws95::iapws95_ideal::*;
use iapws95::iapws95_residual::*;

#[test]
fn debug_sound_speed() {
    let T = 647.0;
    let d = 358.0;
    let w_ref = 252.145078;
    
    let delta = reduced_density(d);
    let tau = inv_reduced_temp(T);
    
    println!("T = {}, d = {}", T, d);
    println!("delta = {}, tau = {}", delta, tau);
    
    let dphi_ddelta = dphi_residual_ddelta(delta, tau);
    let d2phi_ddelta2 = d2phi_residual_ddelta2(delta, tau);
    let d2phi_ddelta_dtau = d2phi_residual_ddelta_dtau(delta, tau);
    let d2phi_dtau2_ideal = d2phi_ideal_dtau2(tau);
    let d2phi_dtau2_residual = d2phi_residual_dtau2(delta, tau);
    
    println!("dphi_ddelta = {}", dphi_ddelta);
    println!("d2phi_ddelta2 = {}", d2phi_ddelta2);
    println!("d2phi_ddelta_dtau = {}", d2phi_ddelta_dtau);
    println!("d2phi_dtau2_ideal = {}", d2phi_dtau2_ideal);
    println!("d2phi_dtau2_residual = {}", d2phi_dtau2_residual);
    
    // Numerical verification of mixed derivative
    let h = 1e-8;
    let dphi_ddelta_tau_plus = dphi_residual_ddelta(delta, tau + h);
    let dphi_ddelta_tau_minus = dphi_residual_ddelta(delta, tau - h);
    let numerical_mixed = (dphi_ddelta_tau_plus - dphi_ddelta_tau_minus) / (2.0 * h);
    println!("Numerical d2phi_ddelta_dtau = {}", numerical_mixed);
    println!("Analytical d2phi_ddelta_dtau = {}", d2phi_ddelta_dtau);
    println!("Difference = {}", (numerical_mixed - d2phi_ddelta_dtau).abs());
    
    // Numerical verification of d2phi_dtau2 (second derivative)
    // f''(τ) ≈ [f'(τ+h) - f'(τ-h)] / (2h) where f' = dphi_residual_dtau
    let dphi_dtau_plus = dphi_residual_dtau(delta, tau + h);
    let dphi_dtau_minus = dphi_residual_dtau(delta, tau - h);
    let numerical_tau2 = (dphi_dtau_plus - dphi_dtau_minus) / (2.0 * h);
    println!("Numerical d2phi_dtau2_residual = {}", numerical_tau2);
    println!("Analytical d2phi_dtau2_residual = {}", d2phi_dtau2_residual);
    println!("Difference = {}", (numerical_tau2 - d2phi_dtau2_residual).abs());
    
    // Numerical verification of d2phi_ddelta2
    let numerical_delta2 = (dphi_residual_ddelta(delta + h, tau) - 2.0 * dphi_residual_ddelta(delta, tau) + dphi_residual_ddelta(delta - h, tau)) / (h * h);
    println!("Numerical d2phi_ddelta2 = {}", numerical_delta2);
    println!("Analytical d2phi_ddelta2 = {}", d2phi_ddelta2);
    println!("Difference = {}", (numerical_delta2 - d2phi_ddelta2).abs());
    
    let numerator = (1.0 + delta * dphi_ddelta - delta * tau * d2phi_ddelta_dtau).powi(2);
    let denominator = tau * tau * (d2phi_dtau2_ideal + d2phi_dtau2_residual);
    
    println!("numerator = {}", numerator);
    println!("denominator = {}", denominator);
    
    let w_squared = IAPWS95_R * T * (
        1.0 + 2.0 * delta * dphi_ddelta + delta * delta * d2phi_ddelta2 
        - numerator / denominator
    );
    
    let w_calc = (w_squared * 1000.0).sqrt();
    
    println!("w_squared = {}", w_squared);
    println!("w_calc = {}", w_calc);
    println!("w_ref = {}", w_ref);
    println!("error = {}", (w_calc - w_ref).abs());
}
