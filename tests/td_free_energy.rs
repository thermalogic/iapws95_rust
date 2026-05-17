#![allow(warnings)]
///  cargo test --test td_free_energy
use assert_approx_eq::assert_approx_eq;
use iapws95::iapws95::*;
use iapws95::iapws95_ideal::*;
use iapws95::iapws95_residual::*;

#[test]
fn test_phi_o_r() {
    // Test Case 1 from Table 6: T = 500 K, rho = 838.025 kg/m³
    let t_test = 500.0;
    let rho_test = 838.025;

    // Reference values from Table 6 (Wagner & Pruss, 2002)
    let phi_o_ref = 2.04797734;     // φ⁰
    let phi_o_d_ref = 0.384236747;  // ∂φ⁰/∂δ
    let phi_o_dd_ref = -0.147637878; // ∂²φ⁰/∂δ²
    let phi_o_t_ref = 9.04611106;   // ∂φ⁰/∂τ
    let phi_o_tt_ref = -1.93249185; // ∂²φ⁰/∂τ²
    let phi_o_dt_ref = 0.0;         // ∂²φ⁰/∂δ∂τ

    let phi_r_ref = -3.42693206;    // φʳ
    let phi_r_d_ref = -0.364366650; // ∂φʳ/∂δ
    let phi_r_dd_ref = 0.856063701; // ∂²φʳ/∂δ²
    let phi_r_t_ref = -5.81403435;  // ∂φʳ/∂τ
    let phi_r_tt_ref = -2.23440737; // ∂²φʳ/∂τ²
    let phi_r_dt_ref = -1.12176915; // ∂²φʳ/∂δ∂τ

    // Calculate our values
    let delta = reduced_density(rho_test);
    let tau = inv_reduced_temp(t_test);

    // Test ideal gas part
    assert_approx_eq!(phi_o_ref, phi_ideal(delta, tau), 1.0e-8f64);
    assert_approx_eq!(phi_o_d_ref, dphi_ideal_ddelta(delta), 1.0e-8f64);
    assert_approx_eq!(phi_o_dd_ref, d2phi_ideal_ddelta2(delta), 1.0e-8f64);
    assert_approx_eq!(phi_o_t_ref, dphi_ideal_dtau(tau), 1.0e-8f64);
    assert_approx_eq!(phi_o_tt_ref, d2phi_ideal_dtau2(tau), 1.0e-8f64);
    assert_approx_eq!(phi_o_dt_ref, d2phi_ideal_dtaudelta(delta, tau), 1.0e-8f64);

    // Test residual part
    assert_approx_eq!(phi_r_ref, phi_residual(delta, tau), 1.0e-8f64);
    assert_approx_eq!(phi_r_d_ref, dphi_residual_ddelta(delta, tau), 1.0e-8f64);
    assert_approx_eq!(phi_r_dd_ref, d2phi_residual_ddelta2(delta, tau), 1.0e-8f64);
    assert_approx_eq!(phi_r_t_ref, dphi_residual_dtau(delta, tau), 1.0e-8f64);
    assert_approx_eq!(phi_r_tt_ref, d2phi_residual_dtau2(delta, tau), 1.0e-8f64);
    assert_approx_eq!(phi_r_dt_ref, d2phi_residual_ddelta_dtau(delta, tau), 1.0e-8f64);
}
