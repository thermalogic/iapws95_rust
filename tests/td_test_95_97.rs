use iapws95::iapws95::*;
use seuif97::{pt, tv, OV,OD, OH, OS, OCP, OCV, OW, OJTC, OIJTC};
use assert_approx_eq::assert_approx_eq;

/// Test structure for thermodynamic property comparison
pub struct TestData {
    pub p: f64,  // Pressure in MPa
    pub t: f64,  // Temperature in c
}

// Test states within IAPWS-95 and IAPWS-IF97 valid range
pub const TEST_STATES: [TestData; 5] = [
    TestData { p: 0.1, t: 300.0 },    // Liquid water
    TestData { p: 0.1, t: 500.0 },    // Superheated steam
    TestData { p: 1.0, t: 500.0 },    // Superheated steam
    TestData { p: 10.0, t: 600.0 },   // High pressure steam
    TestData { p: 20.0, t: 650.0 },   // High pressure, high temperature
];

/// Test structure for new thermodynamic properties
pub struct NewPropData {
    pub t: f64,      // Temperature in Celsius
    pub rho: f64,    // Density in kg/m³
}

// Reference data for testing new properties
// Test states within IAPWS-95 and IAPWS-IF97 valid range
pub const NEW_PROP_DATA: [NewPropData; 1] = [
    NewPropData { t: 300.0-273.15, rho: 0.9965560e3 }, //T: 300.0,    d: 0.9965560e3,   p: 0.992418352e-1,
];

#[test]
fn test_comparison_pt_95_vs_97() {
    // Data flow: IF97 (p,t) -> v and other properties -> rho=1/v -> IAPWS-95 properties -> compare
    println!("\n=== Comparison: IAPWS-95 vs IAPWS-IF97 (p,t domain) ===\n");
    
    for (i, state) in TEST_STATES.iter().enumerate() {
        let p = state.p;
        let t_k = state.t;
        let t_c = t_k - 273.15;
        
        // Step 1: Get specific volume from IF97
        let v_if97 = pt(p, t_c, OV);     
        // Step 2: Calculate density from IF97 specific volume
        let rho = 1.0 / v_if97;
        
        // Step 3: Get properties from IF97
        let h_if97 = pt(p, t_c, OH);
        let s_if97 = pt(p, t_c, OS);
        let cp_if97 = pt(p, t_c, OCP);
        let cv_if97 = pt(p, t_c, OCV);
        let w_if97 = pt(p, t_c, OW);
        let mu_if97 = pt(p, t_c, OJTC);
        
        // Step 4: Get properties from IAPWS-95 using density
        let h_95 = tr2h(t_c, rho);
        let s_95 = tr2s(t_c, rho);
        let cp_95 = tr2cp(t_c, rho);
        let cv_95 = tr2cv(t_c, rho);
        let w_95 = tr2w(t_c, rho);
        //let mu_95 = tr2jt(t_c, rho);
        
              
        // Step 6: Assert approximate equality
        assert_approx_eq!(h_95, h_if97, 5.0e0);
        assert_approx_eq!(s_95, s_if97, 5.0e0);
        assert_approx_eq!(cp_95, cp_if97, 5.0e0);
        assert_approx_eq!(cv_95, cv_if97, 5.0e0);
        assert_approx_eq!(w_95, w_if97, 5.0e0);
        //assert_approx_eq!(mu_95, mu_if97, 1.0);
    }
}

#[test]
fn test_comparison_td_iapws95_vs_if97() {
    // Compare IAPWS-95 results with IAPWS-IF97 (seuif97)
    // Note: These are different formulations, so some differences are expected
    // Data flow: IF97 (t,v) ->  properties，IAPWS-95（t,v) -》properties -> compare
     println!("\n=== Comparison: IAPWS-95 vs IAPWS-IF97 ===\n");
    
    for data in &NEW_PROP_DATA {
        let v = 1.0 / data.rho;
        let t = data.t ;
        // Joule-Thomson coefficient
        let mu_95 = tr2jt(data.t, data.rho);
        let mu_if97 = tv(t, v, OJTC);
        // Isothermal throttling coefficient
        //let mu_t_95 = tr2itt(data.t, data.rho);
        //let mu_t_if97 = tv(t, v, OIJTC);
        assert_approx_eq!( mu_95, mu_if97, 1.0e-1);  
        //assert_approx_eq!( mu_t_95, mu_t_if97, 1.0e-1);  
    }
}


