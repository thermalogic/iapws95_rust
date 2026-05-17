/// Example: Calculate thermodynamic properties of water using IAPWS-95
///
/// This example demonstrates:
/// 1. Single-phase property calculation at given T and rho
/// 2. Saturation property calculation at given T

use iapws95::iapws95::*;
use iapws95::iapws95_saturation::calc_saturation_properties;

fn main() {
    println!("=== IAPWS-95 Thermodynamic Properties Calculator ===\n");

    // Example 1: Single-phase properties at T=500K, rho=838.025 kg/m³
    println!("--- Example 1: Single-phase Properties ---");
    let T = 500.0;
    let rho = 838.025;

    let p = calc_pressure(T, rho);
    let u = calc_internal_energy(T, rho);
    let h = calc_enthalpy(T, rho);
    let s = calc_entropy(T, rho);
    let cv = calc_cv(T, rho);
    let cp = calc_cp(T, rho);
    let w = calc_speed_of_sound(T, rho);

    println!("Input:");
    println!("  Temperature T = {} K", T);
    println!("  Density ρ = {} kg/m³", rho);
    println!("\nOutput:");
    println!("  Pressure p = {:.6} MPa", p);
    println!("  Internal energy u = {:.4} kJ/kg", u);
    println!("  Enthalpy h = {:.4} kJ/kg", h);
    println!("  Entropy s = {:.6} kJ/(kg·K)", s);
    println!("  Cv = {:.6} kJ/(kg·K)", cv);
    println!("  Cp = {:.6} kJ/(kg·K)", cp);
    println!("  Speed of sound w = {:.4} m/s", w);

    // Example 2: Saturation properties at T=450K
    println!("\n--- Example 2: Saturation Properties ---");
    let T_sat = 450.0;

    if let Some(sat) = calc_saturation_properties(T_sat) {
        println!("Input:");
        println!("  Temperature T = {} K", T_sat);
        println!("\nSaturated Liquid ('):");
        println!("  Density ρ' = {:.2} kg/m³", sat.rho_l);
        println!("  Enthalpy h' = {:.2} kJ/kg", sat.h_l);
        println!("  Entropy s' = {:.4} kJ/(kg·K)", sat.s_l);
        println!("\nSaturated Vapor (''):");
        println!("  Density ρ'' = {:.2} kg/m³", sat.rho_v);
        println!("  Enthalpy h'' = {:.2} kJ/kg", sat.h_v);
        println!("  Entropy s'' = {:.4} kJ/(kg·K)", sat.s_v);
        println!("\nSaturation Pressure:");
        println!("  p_sat = {:.6} MPa", sat.p_sat);
    } else {
        println!("Temperature {} K is out of valid saturation range (273.16 - 647.096 K)", T_sat);
    }

    // Example 3: Properties at multiple temperatures
    println!("\n--- Example 3: Saturation Properties at Multiple Temperatures ---");
    println!("{:<10} {:<12} {:<12} {:<12} {:<12} {:<12}", 
             "T(K)", "p_sat(MPa)", "ρ'(kg/m³)", "ρ''(kg/m³)", "h'(kJ/kg)", "h''(kJ/kg)");
    println!("{:-<74}", "");

    for T in [273.16, 373.15, 473.15, 573.15, 647.096] {
        if let Some(sat) = calc_saturation_properties(T) {
            println!("{:<10.2} {:<12.6} {:<12.2} {:<12.2} {:<12.2} {:<12.2}",
                     T, sat.p_sat, sat.rho_l, sat.rho_v, sat.h_l, sat.h_v);
        }
    }
}
