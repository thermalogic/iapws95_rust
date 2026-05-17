/**
 * IAPWS-95 C Interface Example
 * 
 * This example demonstrates how to use the IAPWS-95 library via C bindings.
 * 
 * Compilation:
 *   Windows (MSVC):
 *     cargo build --release
 *     cl /I..\src c_example.c ..\target\release\iapws95.lib
 *   
 *   Linux/Mac:
 *     cargo build --release
 *     gcc -I../src c_example.c -L../target/release -liapws95 -o iapws95_demo
 */

#include <stdio.h>
#include <stdlib.h>
#include "../src/iapws95.h"

int main(void) {
    printf("=== IAPWS-95 Thermodynamic Properties (C Interface) ===\n\n");

    // Get library version
    printf("Library Version: %s\n\n", iapws95_version());

    // Example 2: Saturation properties at T=100°C
    printf("\n--- Example 2: Saturation Properties ---\n");
    double t_sat = 100.0;

    iapws95_saturation_props_t sat;
    int ret = iapws95_saturation_properties(t_sat, &sat);

    if (ret == 0) {
        printf("Input:\n");
        printf("  Temperature T = %.2f °C\n\n", t_sat);
        printf("Saturated Liquid ('):\n");
        printf("  Density ρ' = %.2f kg/m³\n", sat.rho_l);
        printf("  Enthalpy h' = %.2f kJ/kg\n", sat.h_l);
        printf("  Entropy s' = %.4f kJ/(kg·K)\n", sat.s_l);
        printf("\nSaturated Vapor (''):\n");
        printf("  Density ρ'' = %.2f kg/m³\n", sat.rho_v);
        printf("  Enthalpy h'' = %.2f kJ/kg\n", sat.h_v);
        printf("  Entropy s'' = %.4f kJ/(kg·K)\n", sat.s_v);
        printf("\nSaturation Pressure:\n");
        printf("  p_sat = %.6f MPa\n", sat.p_sat);
    } else {
        printf("Error: Temperature out of valid range (0-373.94°C)\n");
    }

    // Example 3: Properties at multiple temperatures
    printf("\n--- Example 3: Saturation Properties Table ---\n");
    printf("%-10s %-12s %-12s %-12s %-12s %-12s\n", 
           "T(°C)", "p_sat(MPa)", "ρ'(kg/m³)", "ρ''(kg/m³)", "h'(kJ/kg)", "h''(kJ/kg)");
    printf("%-74s\n", "--------------------------------------------------------------------------");

    double temps[] = {0.01, 100.0, 200.0, 300.0, 373.94};
    int n_temps = sizeof(temps) / sizeof(temps[0]);

    for (int i = 0; i < n_temps; i++) {
        if (iapws95_saturation_properties(temps[i], &sat) == 0) {
            printf("%-10.2f %-12.6f %-12.2f %-12.2f %-12.2f %-12.2f\n",
                   temps[i], sat.p_sat, sat.rho_l, sat.rho_v, sat.h_l, sat.h_v);
        }
    }

    // Example 5: Direct (T,ρ) → property calculations
    printf("\n--- Example 5: Direct (T,ρ) Property Calculations ---\n");
    double t_direct = 500.0;
    double rho_direct = 838.025;

    double p_calc = iapws95_tr2p(t_direct, rho_direct);
    double u_calc = iapws95_tr2u(t_direct, rho_direct);
    double h_calc = iapws95_tr2h(t_direct, rho_direct);
    double s_calc = iapws95_tr2s(t_direct, rho_direct);
    double cv_calc = iapws95_tr2cv(t_direct, rho_direct);
    double cp_calc = iapws95_tr2cp(t_direct, rho_direct);
    double w_calc = iapws95_tr2w(t_direct, rho_direct);

    printf("Input:\n");
    printf("  Temperature T = %.2f °C\n", t_direct);
    printf("  Density ρ = %.3f kg/m³\n\n", rho_direct);
    printf("Output:\n");
    printf("  Pressure p = %.6f MPa\n", p_calc);
    printf("  Internal energy u = %.4f kJ/kg\n", u_calc);
    printf("  Enthalpy h = %.4f kJ/kg\n", h_calc);
    printf("  Entropy s = %.6f kJ/(kg·K)\n", s_calc);
    printf("  Cv = %.6f kJ/(kg·K)\n", cv_calc);
    printf("  Cp = %.6f kJ/(kg·K)\n", cp_calc);
    printf("  Speed of sound w = %.4f m/s\n", w_calc);

    printf("\n=== Example completed ===\n");

    return 0;
}
