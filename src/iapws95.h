#ifndef IAPWS95_H
#define IAPWS95_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief IAPWS-95 Thermodynamic Properties Library - C Interface
 * 
 * This library provides accurate calculation of thermodynamic properties
 * of water and steam based on the IAPWS-95 formulation.
 * 
 * Valid Range:
 *   Temperature: 273.16 K to 1273 K (0°C to 1000°C)
 *   Pressure: Up to 1000 MPa
 * 
 * Function Categories:
 *   - tr2*: Functions for (t,ρ) → property direct computation
 * 
 * Calling Convention:
 *   All functions use the standard C calling convention (cdecl), which is
 *   compatible with all major C/C++ compilers and platforms. Parameters are
 *   passed on the stack from right to left, and the caller is responsible
 *   for cleaning up the stack. This ensures maximum portability across:
 *   - Windows (MSVC, MinGW)
 *   - Linux (gcc, clang)
 *   - macOS (clang)
 */

 /**
 * @brief Calculate pressure at given temperature and density
 * 
 * @param t_c Temperature in °C
 * @param rho Density in kg/m³
 * @return Pressure p in MPa
 */
double iapws95_tr2p(double t_c, double rho);

/**
 * @brief Calculate specific internal energy at given temperature and density
 * 
 * @param t_c Temperature in °C
 * @param rho Density in kg/m³
 * @return Internal energy u in kJ/kg
 */
double iapws95_tr2u(double t_c, double rho);

/**
 * @brief Calculate specific enthalpy at given temperature and density
 * 
 * @param t_c Temperature in °C
 * @param rho Density in kg/m³
 * @return Enthalpy h in kJ/kg
 */
double iapws95_tr2h(double t_c, double rho);

/**
 * @brief Calculate specific entropy at given temperature and density
 * 
 * @param t_c Temperature in °C
 * @param rho Density in kg/m³
 * @return Entropy s in kJ/(kg·K)
 */
double iapws95_tr2s(double t_c, double rho);

/**
 * @brief Calculate constant-volume specific heat at given temperature and density
 * 
 * @param t_c Temperature in °C
 * @param rho Density in kg/m³
 * @return Cv in kJ/(kg·K)
 */
double iapws95_tr2cv(double t_c, double rho);

/**
 * @brief Calculate constant-pressure specific heat at given temperature and density
 * 
 * @param t_c Temperature in °C
 * @param rho Density in kg/m³
 * @return Cp in kJ/(kg·K)
 */
double iapws95_tr2cp(double t_c, double rho);

/**
 * @brief Calculate speed of sound at given temperature and density
 * 
 * @param t_c Temperature in °C
 * @param rho Density in kg/m³
 * @return Speed of sound w in m/s
 */
double iapws95_tr2w(double t_c, double rho);

// ==========================================================================
// Saturation properties
// ==========================================================================

/**
 * @brief Structure for saturation properties
 */
typedef struct {
    double p_sat;      /**< Saturation vapor pressure (MPa) */
    double rho_l;      /**< Saturated liquid density (kg/m³) */
    double rho_v;      /**< Saturated vapor density (kg/m³) */
    double h_l;        /**< Saturated liquid enthalpy (kJ/kg) */
    double h_v;        /**< Saturated vapor enthalpy (kJ/kg) */
    double s_l;        /**< Saturated liquid entropy (kJ/(kg·K)) */
    double s_v;        /**< Saturated vapor entropy (kJ/(kg·K)) */
} iapws95_saturation_props_t;

/**
 * @brief Calculate all saturation properties at given temperature
 * 
 * @param t_c Temperature in °C
 * @param[out] props Pointer to structure to store results
 * @return 0 on success, -1 on error (temperature out of range)
 */
int iapws95_saturation_properties(double t_c, iapws95_saturation_props_t *props);

/**
 * @brief Get library version string
 * 
 * @return Version string (e.g., "0.1.0")
 */
const char* iapws95_version(void);

#ifdef __cplusplus
}
#endif

#endif /* IAPWS95_H */
