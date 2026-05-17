#![allow(warnings)]
/// cargo test --test T_saturation_table8
/// Test against IAPWS-95 Table 8: Thermodynamic property values in the two-phase region
use iapws95::iapws95_saturation::*;

// Table 8 data from IAPWS-95 (selected temperatures)
// Format: [T(K), p_sat(MPa), rho_l(kg/m³), rho_v(kg/m³), h_l(kJ/kg), h_v(kJ/kg), s_l(kJ/(kg·K)), s_v(kJ/(kg·K))]
const TABLE8_DATA: &[[f64;8]] = &[
    // T,       p_sat,      rho_l,      rho_v,      h_l,        h_v,        s_l,        s_v
    [275.0,     0.698451167e-3,  0.999887406e3,   0.550664919e-2,  0.775972202e1, 0.250428995e4,0.283094670e-1,0.910660121e1],
    [450.0,     0.932203564,    0.890341250e3,   0.481200360e1,   0.749161585e3, 0.277441078e4, 0.210865845e1, 0.660921221e1],
    [625.0,     0.169082693e2,  0.567090385e3,   0.118290280e3,   0.168626976e4, 0.255071625e4,0.380194683e1, 0.518506121e1],
];

fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

#[test]
fn test_table8_saturation_properties() {
    println!("\n{:─^120}", "IAPWS-95 Table 8 Validation");
    println!("{:<8} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}", 
             "T(K)", "p_sat(MPa)", "rho_l(kg/m³)", "rho_v(kg/m³)", "h_l(kJ/kg)", "h_v(kJ/kg)", "s_l(kJ/kg·K)", "s_v(kJ/kg·K)");
    println!("{:─^120}", "");
    
    let mut all_passed = true;
    
    for row in TABLE8_DATA {
        let t = row[0];
        let p_ref = row[1];
        let rho_l_ref = row[2];
        let rho_v_ref = row[3];
        let h_l_ref = row[4];
        let h_v_ref = row[5];
        let s_l_ref = row[6];
        let s_v_ref = row[7];
        
        let sat = calc_saturation_properties(t);
        
        if sat.is_none() {
            println!("T={}: FAILED - calculation returned None", t);
            all_passed = false;
            continue;
        }
        
        let sat = sat.unwrap();
        
        println!("{:<8.1} {:<12.4} {:<12.1} {:<12.2} {:<12.1} {:<12.1} {:<12.4} {:<12.4}", 
                 t, sat.p_sat, sat.rho_l, sat.rho_v, sat.h_l, sat.h_v, sat.s_l, sat.s_v);
        println!("        Ref: {:<12.4} {:<12.1} {:<12.2} {:<12.1} {:<12.1} {:<12.4} {:<12.4}", 
                 p_ref, rho_l_ref, rho_v_ref, h_l_ref, h_v_ref, s_l_ref, s_v_ref);
        
        // Check saturation pressure (within 5% or 0.01 MPa for low pressures)
        let p_tol = if p_ref < 0.01 { 0.001 } else { p_ref * 0.05 };
        if !approx_eq(sat.p_sat, p_ref, p_tol) {
            println!("        ERROR: p_sat mismatch (diff={:.4}, tol={:.4})", (sat.p_sat - p_ref).abs(), p_tol);
            all_passed = false;
        }
        
        // Check liquid density (within 2%)
        if !approx_eq(sat.rho_l, rho_l_ref, rho_l_ref * 0.02) {
            println!("        ERROR: rho_l mismatch (diff={:.1}, tol={:.1})", (sat.rho_l - rho_l_ref).abs(), rho_l_ref * 0.02);
            all_passed = false;
        }
        
        // Check vapor density (within 5%)
        if !approx_eq(sat.rho_v, rho_v_ref, rho_v_ref * 0.05) {
            println!("        ERROR: rho_v mismatch (diff={:.2}, tol={:.2})", (sat.rho_v - rho_v_ref).abs(), rho_v_ref * 0.05);
            all_passed = false;
        }
        
        // Check liquid enthalpy (within 2%)
        if !approx_eq(sat.h_l, h_l_ref, h_l_ref.abs() * 0.02) {
            println!("        ERROR: h_l mismatch (diff={:.1}, tol={:.1})", (sat.h_l - h_l_ref).abs(), h_l_ref.abs() * 0.02);
            all_passed = false;
        }
        
        // Check vapor enthalpy (within 2%)
        if !approx_eq(sat.h_v, h_v_ref, h_v_ref.abs() * 0.02) {
            println!("        ERROR: h_v mismatch (diff={:.1}, tol={:.1})", (sat.h_v - h_v_ref).abs(), h_v_ref.abs() * 0.02);
            all_passed = false;
        }
        
        // Check liquid entropy (within 2%)
        if !approx_eq(sat.s_l, s_l_ref, s_l_ref.abs() * 0.02) {
            println!("        ERROR: s_l mismatch (diff={:.4}, tol={:.4})", (sat.s_l - s_l_ref).abs(), s_l_ref.abs() * 0.02);
            all_passed = false;
        }
        
        // Check vapor entropy (within 2%)
        if !approx_eq(sat.s_v, s_v_ref, s_v_ref.abs() * 0.02) {
            println!("        ERROR: s_v mismatch (diff={:.4}, tol={:.4})", (sat.s_v - s_v_ref).abs(), s_v_ref.abs() * 0.02);
            all_passed = false;
        }
        
        println!();
    }
    
    assert!(all_passed, "Some Table 8 validation tests failed");
}

