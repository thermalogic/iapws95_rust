#![allow(warnings)]
///  cargo test --test td_test
use assert_approx_eq::assert_approx_eq;
use iapws95::iapws95::*;
use iapws95::iapws95_ideal::*;
use iapws95::iapws95_residual::*;

mod td_data;
use td_data::{TD_DATA_TABLE7};

#[test]
fn test_td() {
    for i in 0..11 {
        let T: f64 = TD_DATA_TABLE7[i].T;
        let d: f64 = TD_DATA_TABLE7[i].d;
        let t_c = T - 273.15;
        // Test pressure
        assert_approx_eq!(TD_DATA_TABLE7[i].p, tr2p(t_c, d), 1.0e-6f64);
        assert_approx_eq!(TD_DATA_TABLE7[i].cv, tr2cv(t_c, d), 1.0e-6f64);
        assert_approx_eq!(TD_DATA_TABLE7[i].cv, tr2cv(t_c, d), 1.0e-6f64);
        // Test entropy
        assert_approx_eq!(TD_DATA_TABLE7[i].s, tr2s(t_c, d), 1.0e-6f64);
        // Test speed of sound
        assert_approx_eq!(TD_DATA_TABLE7[i].w, tr2w(t_c, d), 5.0e-6);  
    }
}
