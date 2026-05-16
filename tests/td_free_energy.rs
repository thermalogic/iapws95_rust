#![allow(warnings)]
///  cargo test --test td_free_energy
use assert_approx_eq::assert_approx_eq;
use iapws95::iapws95::*;
use iapws95::iapws95_ideal::*;
use iapws95::iapws95_residual::*;
use std::env;

#[test]
fn test_tv_phi_o_r() {
    // Test case from Table 6: T = 500 K, rho = 838.025 kg/m³
    let t_test = 500.0;
    let rho_test = 838.025;

    // Reference values - from original implementation (n₁ constant, not n₁/τ)
    let phi_o_ref = 2.04797733; // phi0 - ideal gas part
    let phi_r_ref = -3.42693206; // phiR - residual part

    // Calculate our values
    let delta = reduced_density(rho_test);
    let tau = inv_reduced_temp(t_test);

    assert_approx_eq!(phi_o_ref, phi_ideal(delta, tau), 1.0e-6f64);
    assert_approx_eq!(phi_r_ref, phi_residual(delta, tau), 1.0e-6f64);
}
