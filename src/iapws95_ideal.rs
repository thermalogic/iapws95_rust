//! IAPWS-95 Ideal Gas Part Implementation (Internal Module)
//!
//! **This module is internal (`pub(crate)`) and not exposed to external users.**
//!
//! Implements the ideal gas part of the dimensionless Helmholtz free energy
//! based on Equation 5 and Tables 1, 4 of IAPWS-95 Formulation 1995 (Revised 2018).
//!
//! # Formula
//!
//! The ideal gas part of the dimensionless Helmholtz free energy φ°(δ,τ) is given by:
//!
//! ```text
//! φ°(δ,τ) = ln(δ) + n₁° + n₂°τ + n₃°ln(τ) + Σᵢ₌₄⁸ nᵢ°·ln[1 - exp(-γᵢ°τ)]
//! ```
//!
//! where:
//! - δ = ρ/ρc is the reduced density (dimensionless)
//! - τ = Tc/T is the inverse reduced temperature (dimensionless)
//! - nᵢ° and γᵢ° are coefficients from IAPWS-95 Tables 1 and 4
//!
//! # Derivatives
//!
//! The following derivatives are implemented:
//!
//! - First derivative with respect to δ: ∂φ°/∂δ = 1/δ
//! - Second derivative with respect to δ: ∂²φ°/∂δ² = -1/δ²
//! - First derivative with respect to τ: ∂φ°/∂τ = n₂° + n₃°/τ + Σᵢ₌₄⁸ nᵢ°γᵢ°·[(1/[1-exp(-γᵢ°τ)]) - 1]
//! - Second derivative with respect to τ: ∂²φ°/∂τ² = -n₃°/τ² - Σᵢ₌₄⁸ nᵢ°(γᵢ°)²·exp(-γᵢ°τ) / [1-exp(-γᵢ°τ)]²
//! - Mixed derivative: ∂²φ°/∂τ∂δ = 0 (ideal gas has no δ-τ coupling)

// ==========================================================================
// COEFFICIENTS - Ideal Gas Part (Table 1 of IAPWS-95)
// ==========================================================================

/// n coefficients from Table 1 of IAPWS-95 (i=1-8)
const IDEAL_N: [f64; 8] = [
    -8.3204464837497,   // n₁
     6.6832105275932,   // n₂
     3.00632,           // n₃
     0.012436,          // n₄
     0.97315,           // n₅
     1.27950,           // n₆
     0.96956,           // n₇
     0.24873,           // n₈
];

/// γ coefficients for exponential terms (i=4-8) from Table 4 of IAPWS-95
const IDEAL_GAMMA: [f64; 5] = [
    1.28728967,        // γ₄ - for n₄
    3.53734222,        // γ₅ - for n₅
    7.74073708,        // γ₆ - for n₆
    9.24437796,        // γ₇ - for n₇
    27.5075105,        // γ₈ - for n₈
];

// ==========================================================================
// IDEAL-GAS PART CALCULATIONS (Eq. 5 and Table 4 of IAPWS-95)
// ==========================================================================

/// Compute ideal-gas part of dimensionless Helmholtz free energy φ°(δ,τ)
///
/// # Formula
/// ```text
/// φ°(δ,τ) = ln(δ) + n₁ + n₂τ + n₃ln(τ) + Σᵢ₌₄⁸ nᵢ·ln[1 - exp(-γᵢτ)]
/// ```
pub fn phi_ideal(delta: f64, tau: f64) -> f64 {
    let mut sum = delta.ln(); // ln(δ) term

    // Constant and linear terms in τ
    sum += IDEAL_N[0]; // n₁
    sum += IDEAL_N[1] * tau; // n₂τ
    sum += IDEAL_N[2] * tau.ln(); // n₃ln(τ)

    // Exponential terms (i=4-8): Σᵢ₌₄⁸ nᵢln[1-exp(-γᵢτ)]
    // There are 5 exponential terms: n₄,n₅,n₆,n₇,n₈ with γ₄,γ₅,γ₆,γ₇,γ₈
    for i in 0..5 {
        sum += IDEAL_N[3 + i] * (1.0 - (-IDEAL_GAMMA[i] * tau).exp()).ln();
    }

    sum
}

/// Compute first derivative ∂φ°/∂δ = 1/δ
pub fn dphi_ideal_ddelta(delta: f64) -> f64 {
    1.0 / delta
}

/// Compute second derivative ∂²φ°/∂δ² = -1/δ²
pub fn d2phi_ideal_ddelta2(delta: f64) -> f64 {
    -1.0 / (delta * delta)
}

/// Compute first derivative ∂φ°/∂τ = n₂ + n₃/τ + Σᵢ₌₄⁸ nᵢγᵢ·[(1/[1-exp(-γᵢτ)]) - 1]
pub fn dphi_ideal_dtau(tau: f64) -> f64 {
    let mut sum = IDEAL_N[1]; // n₂
    sum += IDEAL_N[2] / tau; // n₃/τ

    // Exponential terms: Σᵢ₌₄⁸ nᵢγᵢ*[(1/[1-exp(-γᵢτ)])-1]
    for i in 0..5 {
        let gamma = IDEAL_GAMMA[i];
        let exp_term = 1.0 / (1.0 - (-gamma * tau).exp());
        sum += IDEAL_N[3 + i] * gamma * (exp_term - 1.0);
    }

    sum
}

/// Compute second derivative ∂²φ°/∂τ² = -n₃/τ² - Σᵢ₌₄⁸ nᵢγᵢ²·exp(-γᵢτ) / [1-exp(-γᵢτ)]²
pub fn d2phi_ideal_dtau2(tau: f64) -> f64 {
    let mut sum = -IDEAL_N[2] / (tau * tau); // -n₃/τ²

    // Exponential terms: -Σᵢ₌₄⁸ nᵢγᵢ²exp(-γᵢτ) / [1-exp(-γᵢτ)]²
    for i in 0..5 {
        let gamma = IDEAL_GAMMA[i];
        let exp_term = (-gamma * tau).exp();
        let denom = 1.0 - exp_term;
        sum -= IDEAL_N[3 + i] * gamma * gamma * exp_term / (denom * denom);
    }

    sum
}

/// Compute mixed derivative ∂²φ°/∂τ∂δ = 0 (ideal gas has no δ-τ coupling)
pub fn d2phi_ideal_dtaudelta(_delta: f64, _tau: f64) -> f64 {
    let _ = _delta;
    let _ = _tau;
    0.0
}
