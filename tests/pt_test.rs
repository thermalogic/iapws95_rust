#![allow(warnings)]
/// cargo test --test pt_test
/// Test inverse solver: (p,T) → ρ using data from td_test.rs
use assert_approx_eq::assert_approx_eq;
use iapws95::iapws95_pT::solve_density;

// Reference data from Table 7 (same as td_test.rs)
// Format: [T(K), rho(kg/m³), p(MPa)]
const PT_data: &[(f64, f64, f64)] = &[
    // T=300K cases (compressed liquid)
    (300.0, 0.9965560e3, 0.992418352e-1),   // p≈0.099 MPa
    (300.0, 0.1005308e4, 0.200022515e2),    // p≈20.0 MPa
    (300.0, 0.1188202e4, 0.700004704e3),    // p≈700 MPa
    // T=500K cases (compressed liquid)
    (500.0, 0.4350000, 0.999679423e-1),     // p≈0.1 MPa (low density vapor)
    (500.0, 0.4532000e1, 0.999938125),      // p≈1.0 MPa (vapor)
    (500.0, 0.8380250e3, 0.100003858e2),    // p≈10 MPa (compressed liquid)
    (500.0, 0.1084564e4, 0.700000405e3),    // p≈700 MPa (compressed liquid)
    // T=647K cases (near critical point)
    (647.0, 0.3580000e3, 0.220384756e2),    // p≈22 MPa (near critical)
    // T=900K cases (high temperature)
    (900.0, 0.2410000, 0.100062559),        // p≈0.1 MPa (vapor)
    (900.0, 0.5261500e2, 0.200000690e2),    // p≈20 MPa (vapor)
    (900.0, 0.8707690e3, 0.700000006e3),    // p≈700 MPa (compressed liquid)
];

#[test]
fn test_pt_to_density() {
    println!("\n{:─^100}", "IAPWS-95 (p,T) → ρ Inverse Solver Validation");
    println!("{:<8} {:<14} {:<14} {:<14} {:<14} {:<12}", 
             "Case", "p(MPa)", "T(K)", "rho_ref(kg/m³)", "rho_calc(kg/m³)", "Error(%)");
    println!("{:-^100}", "");
    
    let mut all_passed = true;
    
    for (i, &(t, p_ref, rho_ref)) in PT_data.iter().enumerate() {
        match solve_density(p_ref, t) {
            Some(rho_calc) => {
                let error_pct = ((rho_calc - rho_ref) / rho_ref * 100.0).abs();
                
                println!("{:<8} {:<14.6} {:<14.1} {:<14.2} {:<14.2} {:<12.4}", 
                         i + 1, p_ref, t, rho_ref, rho_calc, error_pct);
                
                // Check accuracy (within 1% for most cases)
                let tolerance = if t == 647.0 { 0.05 } else { 0.01 };
                if error_pct > tolerance {
                    println!("        ⚠ ERROR: density error {:.4}% exceeds tolerance {:.2}%", 
                             error_pct, tolerance * 100.0);
                    all_passed = false;
                }
            }
            None => {
                println!("{:<8} {:<14.6} {:<14.1} {:<14.2} {:<14} {:<12}", 
                         i + 1, p_ref, t, rho_ref, "None", "-");
                println!("        ✗ FAILED: solve_density returned None");
                all_passed = false;
            }
        }
    }
    
    println!("\n{:─^100}", "");
    assert!(all_passed, "Some (p,T) → ρ inverse solver tests failed");
}

#[test]
fn test_pt_round_trip() {
    println!("\n{:─^80}", "IAPWS-95 Round-Trip Test: ρ → p → ρ");
    println!("{:<8} {:<14} {:<14} {:<14} {:<12}", 
             "Case", "rho_ref", "p_calc", "rho_back", "Error(%)");
    println!("{:-^80}", "");
    
    let mut all_passed = true;
    
    for (i, &(t, rho_ref, p_ref)) in PT_data.iter().enumerate() {
        // Step 1: Calculate pressure from (T, ρ)
        let p_calc = iapws95::iapws95::calc_pressure(t, rho_ref);
        
        // Step 2: Solve density from (p, T)
        match solve_density(p_calc, t) {
            Some(rho_back) => {
                let error_pct = ((rho_back - rho_ref) / rho_ref * 100.0).abs();
                
                println!("{:<8} {:<14.2} {:<14.6} {:<14.2} {:<12.4}", 
                         i + 1, rho_ref, p_calc, rho_back, error_pct);
                
                // Round-trip should be very accurate (< 0.001%)
                let tolerance = 0.001;
                if error_pct > tolerance {
                    println!("        ⚠ ERROR: round-trip error {:.6}% exceeds tolerance {:.3}%", 
                             error_pct, tolerance * 100.0);
                    all_passed = false;
                }
            }
            None => {
                println!("{:<8} {:<14.2} {:<14.6} {:<14} {:<12}", 
                         i + 1, rho_ref, p_calc, "None", "-");
                println!("        ✗ FAILED: round-trip solve_density returned None");
                all_passed = false;
            }
        }
    }
    
    println!("\n{:─^80}", "");
    assert!(all_passed, "Some round-trip tests failed");
}

#[test]
fn test_pt_edge_cases() {
    println!("\n{:─^80}", "IAPWS-95 Edge Cases Test");
    
    // Test invalid inputs
    assert!(solve_density(0.0, 500.0).is_none(), "Should fail for p=0");
    assert!(solve_density(-1.0, 500.0).is_none(), "Should fail for p<0");
    assert!(solve_density(1.0, 0.0).is_none(), "Should fail for T=0");
    assert!(solve_density(1.0, -100.0).is_none(), "Should fail for T<0");
    
    println!("✓ All edge case tests passed (invalid inputs return None)");
}
