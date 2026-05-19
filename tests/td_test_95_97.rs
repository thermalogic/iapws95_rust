#![allow(warnings)]
///  cargo test --test td_test_95_97
use iapws95::iapws95::*;
use seuif97::{pt, tv,OP,OD, OH, OS, OCP, OCV, OW, OJTC, OIJTC,OPC};
use assert_approx_eq::assert_approx_eq;

// IAPWS-IF97
// Table 33. Thermodynamic property values calculated from Eq. (28) for selected values of T and  a
//             T,d,p,h,u,s,cp,w
const r3_Td: [[f64; 8]; 3] = [
    [650., 500., 0.255837018E2, 0.186343019E4, 0.181226279E4, 0.405427273E1, 0.138935717E2, 0.502005554E3],
    [650., 200., 0.222930643E2, 0.237512401E4, 0.226365868E4, 0.485438792E1, 0.446579342E2, 0.383444594E3],
    [750., 500., 0.783095639E2, 0.225868845E4, 0.210206932E4, 0.446971906E1, 0.634165359E1, 0.760696041E3],
];

#[test]
fn test_comparison_Td_95_vs_97() {
    //  IF97 r3 (T,d） compare
    println!("\n=== Comparison: IAPWS-97 vs IAPWS-IF97 (T,t) ===\n");
    
    for i in 0..3 {
        let t_c: f64 = r3_Td[i][0] - 273.15;
        let v: f64 = 1.0 / r3_Td[i][1];
        let rho_95= r3_Td[i][1];
        let p_if97 = tv(t_c, v, OP);
        let h_if97 = tv(t_c, v, OH);
        let s_if97 = tv(t_c,  v, OS);
        let cv_if97 = tv(t_c,  v, OCV);
        let cp_if97 = tv(t_c,  v, OCP);
        let w_if97 = tv(t_c,  v, OW);
        let mu_if97 = tv(t_c, v, OJTC);
        let delta_if97 = tv(t_c, v, OIJTC);
        let beta_if97 = tv(t_c, v, OPC);
        
        //  Get properties from IAPWS-95 using density 
        let p_95 = tr2p(t_c,rho_95);
        let h_95 = tr2h(t_c,rho_95);
        let s_95 = tr2s(t_c, rho_95);
        let cv_95 = tr2cv(t_c, rho_95);
        let cp_95 = tr2cp(t_c, rho_95);
        let w_95 = tr2w(t_c, rho_95);
        let mu_95 = tr2jt(t_c, rho_95);
        let delta_95 = tr2itt(t_c, rho_95);
        let beta_95 = tr2beta_s(t_c, rho_95);
       
        // assert approximate 
        assert_approx_eq!(p_95, p_if97, 1.0e-1);
        assert_approx_eq!(h_95, h_if97, 1.0e0);
        assert_approx_eq!(s_95, s_if97, 1.0e-3);
        assert_approx_eq!(cv_95, cv_if97, 1.0e-1);
        assert_approx_eq!(cp_95, cp_if97, 1.0e0);
        assert_approx_eq!(w_95, w_if97, 1.0e1);
        assert_approx_eq!(mu_95, mu_if97, 1.0e1); 
        assert_approx_eq!(delta_95,delta_if97, 1.0e2);
        assert_approx_eq!(beta_95,beta_if97, 1.0e1);
    }  
 }

