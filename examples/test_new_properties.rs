///  cargo run --example test_new_properties
use iapws95::iapws95::*;

fn main() {
    println!("Testing new thermodynamic property functions\n");
    println!("{:=<80}", "");
    
    // Test cases with different states
    let test_cases = [
        ("Liquid water at 25°C", 25.0, 997.0),
        ("Steam at 100°C (low density)", 100.0, 0.6),
        ("Steam at 200°C", 200.0, 10.0),
        ("Compressed liquid at 350°C", 350.0, 950.0),
        ("High temp steam at 500°C", 500.0, 5.0),
        ("Near critical point", 373.0, 322.0),
    ];
    
    println!("{:<30} {:>12} {:>12} {:>12}", "State", "μ (K/MPa)", "μ_T", "β_s (K/MPa)");
    println!("{:-<80}", "");
    
    for (name, t_c, rho) in test_cases {
        let mu = tr2jt(t_c, rho);
        let mu_T = tr2itt(t_c, rho);
        let beta_s = tr2beta_s(t_c, rho);
        
        println!("{:<30} {:>12.6e} {:>12.6e} {:>12.6e}", name, mu, mu_T, beta_s);
    }
    
    println!("\n{:=<80}", "");
    println!("\nDetailed results:\n");
    
    for (name, t_c, rho) in test_cases {
        let t_k = t_c + 273.15;
        let mu = tr2jt(t_c, rho);
        let mu_T = tr2itt(t_c, rho);
        let beta_s = tr2beta_s(t_c, rho);
        
        println!("State: {}", name);
        println!("  T = {}°C ({} K), ρ = {} kg/m³", t_c, t_k, rho);
        println!("  Joule-Thomson coefficient (μ):          {} K/MPa", mu);
        println!("  Isothermal throttling coefficient (μ_T): {} (dimensionless)", mu_T);
        println!("  Isentropic temp-pressure coeff (β_s):   {} K/MPa", beta_s);
        println!();
    }
}
