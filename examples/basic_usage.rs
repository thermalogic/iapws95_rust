use iapws95::iapws95::*;
use iapws95::iapws95_saturation::sat_t;

fn main() {
    // Single-phase properties at t=300.0-273.15°C, rho= 0.9965560e3 kg/m³
    let t_c = 26.85; // 300.0K;
    let rho = 0.9965560e3;
    println!("Single-phase (t={}°C, ρ={} kg/m³):", t_c, rho);
    println!("  p={:.6} MPa  u={:.4} kJ/kg  h={:.4} kJ/kg  s={:.6} kJ/(kg·K)",
             tr2p(t_c, rho), tr2u(t_c, rho), tr2h(t_c, rho), tr2s(t_c, rho));
    //
    println!("  cv={:.6} kJ/(kg·K)  cp={:.6} kJ/(kg·K)  w={:.4} m/s",
             tr2cv(t_c, rho), tr2cp(t_c, rho), tr2w(t_c, rho));
    //
    println!("  mu={:.6}K/MPa",tr2jt(t_c, rho));
         
    // Saturation properties at t=450.0-273.15°C
    if let Some(sat) = sat_t(450.0-273.15) {
        println!("\nSaturation (t=450.0-273.15°C):");
        println!("  p_sat={:.4} MPa", sat.p_sat);
        println!("  ρ'={:.2} kg/m³  ρ''={:.4} kg/m³", sat.rho_l, sat.rho_v);
        println!("  h'={:.2} kJ/kg  h''={:.2} kJ/kg", sat.h_l, sat.h_v);
        println!("  s'={:.4} kJ/(kg·K)  s''={:.4} kJ/(kg·K)", sat.s_l, sat.s_v);
    }
}
