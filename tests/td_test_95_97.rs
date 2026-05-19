#![allow(warnings)]
///  cargo test --test td_test_95_97
use iapws95::iapws95::*;
use seuif97::{pt, tv,OP,OD, OH, OS, OCP, OCV, OW, OJTC, OIJTC,OPC};
use assert_approx_eq::assert_approx_eq;

mod td_data;
use td_data::{TD_DATA_TABLE7};
#[test]
fn test_comparison_pt_95_vs_97() {
    // Data flow: IF97 (p,t) -> v and other properties -> rho=1/v -> IAPWS-95 properties -> compare
    println!("\n=== Comparison: IAPWS-95 vs IAPWS-IF97 (p,t domain) ===\n");
    
    // Select specific test states: indices 0 (liquid), 3 (superheated steam), 7 (near critical)
    let selected_indices = [3];
    
    for &i in &selected_indices {
        let state = &TD_DATA_TABLE7[i];
        let p_95 = state.p;
        let t_k = state.T;
        let rho_95 = state.d;
        let v_95=1.0/rho_95;
        let t_c = t_k - 273.15;
       
        // Step 1: Getproperties from IF97
        let rho_if97 = pt(p_95, t_c, OD);     
        let p_if97 = tv(t_c, v_95, OP);
        let h_if97 = tv(t_c, v_95, OH);
        let s_if97 = tv(t_c,  v_95, OS);
        let cv_if97 = tv(t_c,  v_95, OCV);
        let cp_if97 = tv(t_c,  v_95, OCP);
        let w_if97 = tv(t_c,  v_95, OW);
        let mu_if97 = tv(t_c, v_95, OJTC);
        let delta_if97 = tv(t_c, v_95, OIJTC);
        let beta_if97 = tv(t_c, v_95, OPC);
        
        // Step 2: Get properties from IAPWS-95 using density 
        let h_95 = tr2h(t_c,rho_95);
        let s_95 = tr2s(t_c, rho_95);
        let cv_95 = tr2cv(t_c, rho_95);
        let cp_95 = tr2cp(t_c, rho_95);
        let w_95 = tr2w(t_c, rho_95);
        let mu_95 = tr2jt(t_c, rho_95);
        let delta_95 = tr2itt(t_c, rho_95);
        let beta_95 = tr2beta_s(t_c, rho_95);
       
        // Step 3 Assert approximate 
        assert_approx_eq!(p_95, p_if97, 5.0e-2);
        assert_approx_eq!(rho_95, rho_if97, 5.0e1);
        assert_approx_eq!(h_95, h_if97, 5.0e-1);
        assert_approx_eq!(s_95, s_if97, 5.0e-1);
        assert_approx_eq!(cv_95, cv_if97, 5.0e-1);
        assert_approx_eq!(cp_95, cp_if97, 5.0e-1);
        assert_approx_eq!(w_95, w_if97, 5.0e-1);
        // TODO: left: `0.016006260491478408`, right: `19.138216622369853
        assert_approx_eq!(mu_95, mu_if97, 1.0e-2); 
        assert_approx_eq!(delta_95,delta_if97, 1.0e1);
        assert_approx_eq!(beta_95,beta_if97, 1.0e1);
    }
}


