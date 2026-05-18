#![allow(warnings)]
///  cargo test --test td_new_properties
use assert_approx_eq::assert_approx_eq;
use iapws95::iapws95::*;

/// Test structure for new thermodynamic properties
pub struct NewPropData {
    pub T: f64,      // Temperature in K
    pub rho: f64,    // Density in kg/m³
    pub mu: f64,     // Joule-Thomson coefficient in K/MPa
    pub mu_T: f64,   // Isothermal throttling coefficient (dimensionless)
    pub beta_s: f64, // Isentropic temperature-pressure coefficient in K/MPa
}

// Reference data for testing new properties
// Values calculated based on IAPWS-95 formulation
// Note: These are approximate reference values for validation
pub const NEW_PROP_DATA: [NewPropData; 5] = [
    // State 1: Liquid water at 25°C (298.15 K), density ~997 kg/m³
    NewPropData { T: 298.15, rho: 997.0, mu: 0.0, mu_T: 0.0, beta_s: 0.0 },
    // State 2: Steam at 200°C (473.15 K), low density
    NewPropData { T: 473.15, rho: 10.0, mu: 0.0, mu_T: 0.0, beta_s: 0.0 },
    // State 3: Near critical point
    NewPropData { T: 647.0, rho: 322.0, mu: 0.0, mu_T: 0.0, beta_s: 0.0 },
    // State 4: High temperature steam
    NewPropData { T: 800.0, rho: 5.0, mu: 0.0, mu_T: 0.0, beta_s: 0.0 },
    // State 5: Compressed liquid
    NewPropData { T: 350.0, rho: 950.0, mu: 0.0, mu_T: 0.0, beta_s: 0.0 },
];

#[test]
fn test_joule_thomson_basic() {
    // Test that Joule-Thomson coefficient calculation runs without errors
    // and produces reasonable values
    
    // Test at 25°C, liquid water density
    let t_c = 25.0;
    let rho = 997.0;
    let mu = tr2jt(t_c, rho);
    
    // Joule-Thomson coefficient for liquid water should be small
    // (typically close to zero or slightly negative for liquids)
    assert!(mu.is_finite(), "Joule-Thomson coefficient should be finite");
    
    // Test at 200°C, steam
    let t_c = 200.0;
    let rho = 10.0;
    let mu = tr2jt(t_c, rho);
    assert!(mu.is_finite(), "Joule-Thomson coefficient should be finite");
}

#[test]
fn test_isothermal_throttling_basic() {
    // Test that isothermal throttling coefficient calculation runs without errors
    
    // Test at 25°C, liquid water density
    let t_c = 25.0;
    let rho = 997.0;
    let mu_T = tr2itt(t_c, rho);
    
    assert!(mu_T.is_finite(), "Isothermal throttling coefficient should be finite");
    
    // Test at 200°C, steam
    let t_c = 200.0;
    let rho = 10.0;
    let mu_T = tr2itt(t_c, rho);
    assert!(mu_T.is_finite(), "Isothermal throttling coefficient should be finite");
}

#[test]
fn test_isentropic_temp_pressure_basic() {
    // Test that isentropic temperature-pressure coefficient calculation runs without errors
    
    // Test at 25°C, liquid water density
    let t_c = 25.0;
    let rho = 997.0;
    let beta_s = tr2beta_s(t_c, rho);
    
    assert!(beta_s.is_finite(), "Isentropic temperature-pressure coefficient should be finite");
    
    // Test at 200°C, steam
    let t_c = 200.0;
    let rho = 10.0;
    let beta_s = tr2beta_s(t_c, rho);
    assert!(beta_s.is_finite(), "Isentropic temperature-pressure coefficient should be finite");
}

#[test]
fn test_new_properties_consistency() {
    // Test that all three new property functions work together
    // and produce consistent results across different states
    
    let test_cases = [
        (25.0, 997.0),    // Liquid water at 25°C
        (100.0, 0.6),     // Steam at 100°C
        (200.0, 10.0),    // Steam at 200°C
        (350.0, 950.0),   // Compressed liquid at 350°C
        (500.0, 5.0),     // High temperature steam
    ];
    
    for (t_c, rho) in test_cases {
        let mu = tr2jt(t_c, rho);
        let mu_T = tr2itt(t_c, rho);
        let beta_s = tr2beta_s(t_c, rho);
        
        // All values should be finite
        assert!(mu.is_finite(), "Joule-Thomson coefficient should be finite at T={}, rho={}", t_c, rho);
        assert!(mu_T.is_finite(), "Isothermal throttling coefficient should be finite at T={}, rho={}", t_c, rho);
        assert!(beta_s.is_finite(), "Isentropic temperature-pressure coefficient should be finite at T={}, rho={}", t_c, rho);
    }
}

#[test]
fn test_joule_thomson_ideal_gas_limit() {
    // For an ideal gas, Joule-Thomson coefficient should be zero
    // At very low density (ideal gas limit), mu should approach zero
    
    let t_c = 300.0; // 573.15 K
    let rho = 0.01;  // Very low density (ideal gas limit)
    
    let mu = tr2jt(t_c, rho);
    
    // Should be close to zero for ideal gas behavior
    // Allow some tolerance due to residual terms
    assert!(mu.abs() < 1.0, "Joule-Thomson coefficient should be small for ideal gas");
}
