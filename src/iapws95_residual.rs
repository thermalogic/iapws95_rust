//! IAPWS-95 Residual Part Implementation
//!
//! Implements the residual part of the dimensionless Helmholtz free energy
//! based on Table 5 of IAPWS-95 Formulation 1995 (Revised 2018).
//!
//! # Formula
//!
//! ```text
//! φʳ(δ,τ) = Σᵢ nᵢ·δᵈⁱ·τᵗⁱ                                    [polynomial terms, i=1-7]
//!         + Σᵢ nᵢ·δᵈⁱ·τᵗⁱ·exp(-δᶜⁱ)                          [exponential terms, i=8-51]
//!         + Σᵢ nᵢ·δᵈⁱ·τᵗⁱ·exp[-αᵢ(δ-εᵢ)²-βᵢ(τ-γᵢ)²]          [Gaussian terms, i=52-54]
//!         + Σᵢ nᵢ·Δᵇⁱ·δ·F(δ,τ)                               [non-analytic terms, i=55-56]
//! ```

// ==========================================================================
// COEFFICIENTS - Residual Part (Table 5 of IAPWS-95)
// ==========================================================================

/// Polynomial term: (d exponent of delta, t exponent of tau, n coefficient)
struct PolyTerm {
    d: i32,
    t: f64,
    n: f64,
}

/// Exponential term with c=1 or c=2: (d exponent of delta, t exponent of tau, n coefficient)
struct ExpTermC1C2 {
    d: i32,
    t: i32,
    n: f64,
}

/// Exponential term with variable c: (c exponent of delta, d exponent of delta, t exponent of tau, n coefficient)
struct ExpTermCN {
    c: i32,
    d: i32,
    t: i32,
    n: f64,
}

/// Gaussian term: (d, t, n, alpha, beta, gamma, epsilon)
struct GaussTerm {
    d: i32,
    t: i32,
    n: f64,
    a: i32,  // alpha
    b: i32,  // beta
    g: f64,  // gamma
    e: i32,  // epsilon
}

/// Non-analytic term: (alpha, beta, B, n, C, D, A, bt)
struct NonAnalTerm {
    a: f64,  // alpha for Delta calculation
    b: f64,  // beta exponent
    B: f64,
    n: f64,
    C: i32,
    D: i32,
    A: f64,
    bt: f64,
}

// Polynomial terms (i=1 to i=7): nᵢδᵈⁱτᵗⁱ
const RES_POLY_D1: [PolyTerm; 7] = [
    PolyTerm { d: 1, t: -0.5, n: 0.12533547935523e-1 },   // i=1
    PolyTerm { d: 1, t:  0.875, n: 0.78957634722828e1 },   // i=2
    PolyTerm { d: 1, t:  1.0, n: -0.87803203303561e1 },   // i=3
    PolyTerm { d: 2, t:  0.5, n: 0.31802509345418 },     // i=4
    PolyTerm { d: 2, t:  0.75, n: -0.26145533859358 },     // i=5
    PolyTerm { d: 3, t:  0.375, n: -0.78199751687981e-2 },  // i=6
    PolyTerm { d: 4, t:  1.0, n: 0.88089493102134e-2 },  // i=7
];

// Exponential terms (i=8 to i=22): nᵢδᵈⁱτᵗⁱexp(-δ) with c=1 for all
const RES_EXP_D2_C1: [ExpTermC1C2; 15] = [
    ExpTermC1C2 { d: 1, t:   4, n: -0.66856572307965 },       // i=1
    ExpTermC1C2 { d: 1, t:   6, n:  0.20433810950965 },       // i=2
    ExpTermC1C2 { d: 1, t:  12, n: -0.66212605039687e-4 },    // i=3
    ExpTermC1C2 { d: 2, t:   1, n: -0.19232721156002 },       // i=4
    ExpTermC1C2 { d: 2, t:   5, n: -0.25709043003438 },       // i=5
    ExpTermC1C2 { d: 3, t:   4, n:  0.16074868486251 },       // i=6
    ExpTermC1C2 { d: 4, t:   2, n: -0.40092828925807e-1 },    // i=7
    ExpTermC1C2 { d: 4, t:  13, n:  0.39343422603254e-6 },    // i=8
    ExpTermC1C2 { d: 5, t:   9, n: -0.75941377088144e-5 },    // i=9
    ExpTermC1C2 { d: 7, t:   3, n:  0.56250979351888e-3 },    // i=10
    ExpTermC1C2 { d: 9, t:   4, n: -0.15608652257135e-4 },    // i=11
    ExpTermC1C2 { d: 10, t: 11, n:  0.11537996422951e-8 },    // i=12
    ExpTermC1C2 { d: 11, t:  4, n:  0.36582165144204e-6 },    // i=13
    ExpTermC1C2 { d: 13, t: 13, n: -0.13251180074668e-11 },   // i=14
    ExpTermC1C2 { d: 15, t:  1, n: -0.62639586912454e-9 },    // i=15
];

// Exponential terms (i=23 to i=42): nᵢδᵈⁱτᵗⁱexp(-δ*δ) with c=2 for all
const RES_EXP_D2_C2: [ExpTermC1C2; 20] = [
    ExpTermC1C2 { d: 1, t:   7, n: -0.10793600908932 },       // i=16
    ExpTermC1C2 { d: 2, t:   1, n:  0.17611491008752e-1 },    // i=17
    ExpTermC1C2 { d: 2, t:   9, n:  0.22132295167546 },       // i=18
    ExpTermC1C2 { d: 2, t:  10, n: -0.40247669763528 },       // i=19
    ExpTermC1C2 { d: 3, t:  10, n:  0.58083399985759 },       // i=20
    ExpTermC1C2 { d: 4, t:   3, n:  0.49969146990806e-2 },    // i=21
    ExpTermC1C2 { d: 4, t:   7, n: -0.31358700712549e-1 },    // i=22
    ExpTermC1C2 { d: 4, t:  10, n: -0.74315929710341 },       // i=23
    ExpTermC1C2 { d: 5, t:  10, n:  0.47807329915480 },       // i=24
    ExpTermC1C2 { d: 6, t:   6, n:  0.20527940895948e-1 },    // i=25
    ExpTermC1C2 { d: 6, t:  10, n: -0.13636435110343 },       // i=26
    ExpTermC1C2 { d: 7, t:  10, n:  0.14180634400617e-1 },    // i=27
    ExpTermC1C2 { d: 9, t:   1, n:  0.83326504880713e-2 },    // i=28
    ExpTermC1C2 { d: 9, t:   2, n: -0.29052336009585e-1 },    // i=29
    ExpTermC1C2 { d: 9, t:   3, n:  0.38615085574206e-1 },    // i=30
    ExpTermC1C2 { d: 9, t:   4, n: -0.20393486513704e-1 },    // i=31
    ExpTermC1C2 { d: 9, t:   8, n: -0.16554050063734e-2 },    // i=32
    ExpTermC1C2 { d: 10, t:  6, n:  0.19955571979541e-2 },    // i=33
    ExpTermC1C2 { d: 10, t:  9, n:  0.15870308324157e-3 },    // i=34
    ExpTermC1C2 { d: 12, t:  8, n: -0.16388568342530e-4 },    // i=35
];

// Exponential terms (i=43 to i=51): nᵢδᵈⁱτᵗⁱexp(-γᵢδᶜⁱ) with c=3,4,6 for all
const RES_EXP_D2_CN: [ExpTermCN; 9] = [
    // i=43 to i=46 (c=3)
    ExpTermCN { c: 3, d: 3, t: 16, n:  0.43613615723811e-1 },    // i=43
    ExpTermCN { c: 3, d: 4, t: 22, n:  0.34994005463765e-1 },    // i=44
    ExpTermCN { c: 3, d: 4, t: 23, n: -0.76788197844621e-1 },    // i=45
    ExpTermCN { c: 3, d: 5, t: 23, n:  0.22446277332006e-1 },    // i=46

    // i=47 (c=4)
    ExpTermCN { c: 4, d: 14, t: 10, n: -0.62689710414685e-4 },   // i=47

    // i=48 to i=51 (c=6)
    ExpTermCN { c: 6, d: 3, t: 50, n: -0.55711118565645e-9 },    // i=48
    ExpTermCN { c: 6, d: 6, t: 44, n: -0.19905718354408 },       // i=49
    ExpTermCN { c: 6, d: 6, t: 46, n:  0.31777497330738 },       // i=50
    ExpTermCN { c: 6, d: 6, t: 50, n: -0.11841182425981 },       // i=51
];

/// Gaussian terms (i=52-54) from Table 5 of IAPWS-95
const RES_GAUSS: [GaussTerm; 3] = [
    GaussTerm { d: 3, t: 0, n: -0.31306260323435e2, a: 20, b: 150, g: 1.21, e: 1 },   // i=1
    GaussTerm { d: 3, t: 1, n:  0.31546140237781e2, a: 20, b: 150, g: 1.21, e: 1 },   // i=2
    GaussTerm { d: 3, t: 4, n: -0.25213154341695e4, a: 20, b: 250, g: 1.25, e: 1 },   // i=3
];

/// Non-analytic terms (i=55-56) from Table 5 of IAPWS-95: nᵢ·Δᵇⁱ·δ·F(δ,τ)
const RES_NON_ANAL: [NonAnalTerm; 2] = [
    NonAnalTerm { a: 3.5, b: 0.85, B: 0.2, n: -0.14874640856724, C: 28, D: 700, A: 0.32, bt: 0.3 },   // i=1
    NonAnalTerm { a: 3.5, b: 0.95, B: 0.2, n:  0.31806110878444, C: 32, D: 800, A: 0.32, bt: 0.3 },   // i=2
];

// ==========================================================================
// HELPER - Gaussian term exponential argument
// ==========================================================================

/// Compute exp[-α(δ-ε)²-β(τ-γ)²] for Gaussian term
#[inline]
fn gauss_exp(term: &GaussTerm, delta: f64, tau: f64) -> f64 {
    let d_e = delta - term.e as f64;
    let t_g = tau - term.g;
    (-(term.a as f64) * d_e * d_e - (term.b as f64) * t_g * t_g).exp()
}

// ==========================================================================
// RESIDUAL PART CALCULATIONS (Eq. 6 and Table 5 of IAPWS-95)
// ==========================================================================

/// Compute residual part of dimensionless Helmholtz free energy φʳ(δ,τ)
///
/// Sums all 56 terms across four categories: polynomial (7), exponential (44),
/// Gaussian (3), and non-analytic (2).
pub fn phi_residual(delta: f64, tau: f64) -> f64 {
    let mut sum = 0.0f64;

    // Polynomial terms (i=1 to i=7): Σᵢ nᵢδᵈⁱτᵗⁱ
    for term in &RES_POLY_D1 {
        sum += term.n * delta.powi(term.d) * tau.powf(term.t);
    }

    // Exponential terms (i=8 to 22): Σᵢ nᵢδᵈⁱτᵗⁱexp(-δ) c=1
    for term in &RES_EXP_D2_C1 {
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * (-delta).exp();
    }

    // Exponential terms (i=23 to 42): Σᵢ nᵢδᵈⁱτᵗⁱexp(-δ*δ) c=2
    let delta_2 = delta * delta;
    for term in &RES_EXP_D2_C2 {
        sum += term.n *delta.powi(term.d) * tau.powi(term.t) * (-delta_2).exp();
    }

    // Exponential terms (i=43 to 51): Σᵢ nᵢδᵈⁱτᵗⁱexp(-δᶜⁱ) c=3,4,6
    let delta_3 = delta_2 * delta;
    let delta_4 = delta_3 * delta;
    let delta_6 = delta_3 * delta_3;

    // c=3 terms (i=43 to i=46)
    for term in &RES_EXP_D2_CN[0..4] {
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * (-delta_3).exp();
    }

    // c=4 term (i=47)
    let term = &RES_EXP_D2_CN[4];
    sum += term.n * delta.powi(term.d) * tau.powi(term.t) * (-delta_4).exp();

    // c=6 terms (i=48 to i=51)
    for term in &RES_EXP_D2_CN[5..9] {
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * (-delta_6).exp();
    }

    // Gaussian terms (i=52 to i=54): Σᵢ nᵢδᵈⁱτᵗⁱexp[-αᵢ(δ-εᵢ)²-βᵢ(τ-γᵢ)²]
    for term in &RES_GAUSS {
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * gauss_exp(term, delta, tau);
    }

    // Non-analytic terms (i=55 to i=56): Σᵢ nᵢΔᵇⁱδF(δ,τ)
    for term in &RES_NON_ANAL {
        let d_1 = delta - 1.0;
        let d_1_2 = d_1 * d_1;
        let tita = (1.0 - tau) + term.A * d_1_2.powf(0.5 / term.bt);
        let term2 = tau - 1.0;
        let f_val = (-term.C as f64 * d_1_2 - term.D as f64 * term2 * term2).exp();
        let delta_val = tita * tita + term.B * d_1_2.powf(term.a);
        sum += term.n * delta_val.powf(term.b) * delta * f_val;
    }

    sum
}

/// Compute first derivative ∂φʳ/∂δ for residual part
///
/// Applies analytical derivatives to each term category. For non-analytic terms,
/// applies the chain rule to compute both direct δ derivative and indirect contribution through Δ(δ,τ).
pub fn dphi_residual_ddelta(delta: f64, tau: f64) -> f64 {
    let mut sum = 0.0f64;

    // Polynomial terms (i=1 to i=7): ∂(Σᵢ nᵢδᵈⁱτᵗⁱ)/∂δ = nᵢdᵢδᵈⁱ⁻¹τᵗⁱ
    for term in &RES_POLY_D1 {
        sum += term.n * (term.d as f64) * delta.powi(term.d - 1) * tau.powf(term.t);
    }

    // Exponential terms (i=8 to i=22): c=1 ∂(Σᵢ nᵢδᵈⁱτᵗⁱexp(-δᶜⁱ))/∂δ
    for term in &RES_EXP_D2_C1 {
        let delta_c = delta;
        let exp_term = (-delta_c).exp();
        let deriv = (term.d as f64) - delta_c;
        sum += term.n * exp_term * delta.powi(term.d - 1) *tau.powi(term.t) * deriv;
    }

    // Exponential terms (i=23 to i=42): c=2 ∂(Σᵢ nᵢδᵈⁱτᵗⁱexp(-δ²))/∂δ
    for term in &RES_EXP_D2_C2 {
        let delta_c = delta * delta;
        let exp_term = (-delta_c).exp();
        let deriv = (term.d as f64) - 2.0 * delta_c;
        sum += term.n * exp_term * delta.powi(term.d - 1) * tau.powi(term.t) * deriv;
    }

    // Exponential terms (i=43 to i=51): c=3,4,6 ∂(Σᵢ nᵢδᵈⁱτᵗⁱexp(-δᶜⁱ))/∂δ
    for term in &RES_EXP_D2_CN {
        let delta_c = delta.powi(term.c);
        let exp_term = (-delta_c).exp();
        let deriv = (term.d as f64) - (term.c as f64) * delta_c;
        sum += term.n * exp_term * delta.powi(term.d - 1) * tau.powi(term.t) * deriv;
    }

    // Gaussian terms (i=1 to i=3): ∂(Σᵢ nᵢδᵈⁱτᵗⁱexp[-αᵢ(δ-εᵢ)²-βᵢ(τ-γᵢ)²])/∂δ
    for term in &RES_GAUSS {
        let exp_term = gauss_exp(term, delta, tau);
        let deriv = (term.d as f64) / delta - 2.0 * (term.a as f64) * (delta - term.e as f64);
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * exp_term * deriv;
    }

    // Non-analytic terms (i=1 to i=2): ∂(Σᵢ nᵢΔᵇⁱδF)/∂δ
    for term in &RES_NON_ANAL {
        let d_1 = delta - 1.0;
        let d_1_2 = d_1 * d_1;
        let tita = (1.0 - tau) + term.A * d_1_2.powf(0.5 / term.bt);
        let f_val = (-term.C as f64 * d_1_2 - term.D as f64 * (tau - 1.0).powi(2)).exp();
        let f_d = -2.0 * term.C as f64 * f_val * d_1;

        let tita2 = tita * tita;
        let delta_val = tita2 + term.B * d_1_2.powf(term.a);
        let delta_d = d_1
            * (term.A * tita * 2.0 / term.bt * d_1_2.powf(0.5 / term.bt - 1.0)
                + 2.0 * term.B * term.a * d_1_2.powf(term.a - 1.0));

        let (delta_b, delta_bd) = if delta_val == 0.0 {
            (0.0, 0.0)
        } else {
            let db = delta_val.powf(term.b);
            let dbd = term.b * delta_val.powf(term.b - 1.0) * delta_d;
            (db, dbd)
        };

        sum += term.n * (delta_b * (f_val + delta * f_d) + delta_bd * delta * f_val);
    }

    sum
}

/// Compute second derivative ∂²φʳ/∂δ² for residual part
pub fn d2phi_residual_ddelta2(delta: f64, tau: f64) -> f64 {
    let mut sum = 0.0f64;

    // Polynomial terms (i=1 to i=7): ∂²(nᵢ·δᵈⁱ·τᵗⁱ)/∂δ² = nᵢ·dᵢ·(dᵢ-1)·δᵈⁱ⁻²·τᵗⁱ
    for term in &RES_POLY_D1 {
        let d = term.d as f64;
        if term.d >= 2 {
            let contrib = term.n * d * (d - 1.0) * delta.powi(term.d - 2) * tau.powf(term.t);
            sum += contrib;
        }
    }

    // Exponential terms (i=8 to i=22): nᵢ δᵈⁱ τᵗⁱ exp(-δ), c= 1,fn error
    for term in &RES_EXP_D2_C1 {
        let d = term.d as f64;
        let exp_term: f64 = (-delta).exp();
        let d_term = delta.powi(term.d-2)*tau.powi(term.t)*((d - delta)*(d - 1.0-delta)-delta);
        sum += term.n*exp_term*d_term;
    }

    // Exponential terms (i=23 to i=42): nᵢ exp(-δ²)[δᵈⁱ τᵗⁱ()] , c=2 fn error
    for term in &RES_EXP_D2_C2 {
        let d = term.d as f64;
        let delta_c= delta *delta;
        let exp_term: f64 = (-delta_c).exp();
        let d_term = delta.powi(term.d-2)*tau.powi(term.t)*((d -2.0*delta_c)*(d - 1.0-delta_c)-2.0*delta_c);
        sum += term.n*exp_term*d_term;
       }

    // Exponential terms (i=43 to i=51): nᵢ  exp(-δᶜⁱ)[δᵈⁱ τᵗⁱ()] c=3,4,6 fn error
    for term in &RES_EXP_D2_CN {
        let d = term.d as f64;
        let delta_c = delta.powi(term.c);
        let exp_term: f64 = (-delta_c).exp();
        let d_term = delta.powi(term.d-2)*tau.powi(term.t)*((d -2.0*delta_c)*(d - 1.0-delta_c)-2.0*delta_c);
        sum += term.n*exp_term*d_term;      
    }

    // Gaussian terms
    for term in &RES_GAUSS {
        let d = term.d as f64;
        let a = term.a as f64;
        let e = term.e as f64;
        
        let delta_e = delta - e;
        let exp_term = gauss_exp(term, delta, tau);
        
        // Second derivative of δᵗ * exp[-α(δ-ε)²]
        // g(δ) = δᵈ * exp[-α(δ-ε)²]
        // g'(δ) = δᵈ * exp[-α(δ-ε)²] * (d/δ - 2α(δ-ε))
        // g''(δ) = δᵈ * exp[-α(δ-ε)²] * [(d/δ - 2α(δ-ε))² + (-d/δ² - 2α)]
        let first_deriv_factor = d / delta - 2.0 * a * delta_e;
        let second_deriv_correction = -d / (delta * delta) - 2.0 * a;
        let factor = first_deriv_factor * first_deriv_factor + second_deriv_correction;
        
        let contrib = term.n * delta.powi(term.d) * tau.powi(term.t) * exp_term * factor;
        sum += contrib;
    }

    // Now, let's add Non-analytic terms!
    for term in &RES_NON_ANAL {
        let d_1 = delta - 1.0;
        let d_1_2 = d_1 * d_1;
        let tita = (1.0 - tau) + term.A * d_1_2.powf(0.5 / term.bt);
        let f_val = (-term.C as f64 * d_1_2 - term.D as f64 * (tau - 1.0).powi(2)).exp();
        
        let f_d = -2.0 * term.C as f64 * f_val * d_1;
        let f_dd = -2.0 * term.C as f64 * (f_val + d_1 * f_d);

        let tita2 = tita * tita;
        let delta_val = tita2 + term.B * d_1_2.powf(term.a);
        
        let delta_d = d_1 * (
            term.A * tita * 2.0 / term.bt * d_1_2.powf(0.5 / term.bt - 1.0) 
            + 2.0 * term.B * term.a * d_1_2.powf(term.a - 1.0)
        );
        
        // Compute delta_dd (second derivative of delta_val with respect to delta)
        let tita_d = term.A * d_1 * 2.0 / term.bt * d_1_2.powf(0.5 / term.bt - 1.0);
        let tita_dd = term.A * (2.0 / term.bt) * (
            d_1_2.powf(0.5 / term.bt - 1.0) 
            + d_1 * (0.5 / term.bt - 1.0) * d_1_2.powf(0.5 / term.bt - 2.0) * 2.0 * d_1
        );
        let delta_dd = 2.0 * tita_d * tita_d + 2.0 * tita * tita_dd 
            + term.B * (
                2.0 * term.a * d_1_2.powf(term.a - 1.0) 
                + term.a * (term.a - 1.0) * d_1_2.powf(term.a - 2.0) * 4.0 * d_1_2
            );

        let (delta_b, delta_bd, delta_bdd) = if delta_val == 0.0 {
            (0.0, 0.0, 0.0)
        } else {
            let db = delta_val.powf(term.b);
            let dbd = term.b * delta_val.powf(term.b - 1.0) * delta_d;
            let dbdd = term.b * (
                (term.b - 1.0) * delta_val.powf(term.b - 2.0) * delta_d * delta_d
                + delta_val.powf(term.b - 1.0) * delta_dd
            );
            (db, dbd, dbdd)
        };

        // Now compute second derivative of the term nᵢ * (Δ^bᵢ * δ * F)
        let term_dd = term.n * (
            delta_bdd * delta * f_val
            + 2.0 * delta_bd * (f_val + delta * f_d)
            + delta_b * (2.0 * f_d + delta * f_dd)
        );
        sum += term_dd;
    }

    sum
}

/// Compute first derivative ∂φʳ/∂τ for residual part
pub fn dphi_residual_dtau(delta: f64, tau: f64) -> f64 {
    let mut sum = 0.0f64;

    // Polynomial terms (i=1 to i=7): ∂(Σᵢ nᵢ·δᵈⁱ·τᵗⁱ)/∂τ = nᵢ·tᵢ·δᵈⁱ·τᵗⁱ⁻¹
    for term in &RES_POLY_D1 {
        sum += term.n * term.t * delta.powi(term.d) * tau.powf(term.t - 1.0);
    }

    // Exponential terms (i=8 to i=22): c=1 ∂(Σᵢ nᵢδᵈⁱτᵗⁱexp(-δ))/∂τ
    for term in &RES_EXP_D2_C1 {
        sum += term.n * (term.t as f64) * delta.powi(term.d) * tau.powf((term.t - 1) as f64)
            * (-delta).exp();
    }

    // Exponential terms (i=23 to i=42): c=2 ∂(Σᵢ nᵢδᵈⁱτᵗⁱexp(-δ²))/∂τ
    let delta_2 = delta * delta;
    for term in &RES_EXP_D2_C2 {
        sum += term.n * (term.t as f64) * delta.powi(term.d) * tau.powf((term.t - 1) as f64)
            * (-delta_2).exp();
    }

    // Exponential terms (i=43 to i=51): c=3,4,6 ∂(Σᵢ nᵢδᵈⁱτᵗⁱexp(-δᶜ))/∂τ
    for term in &RES_EXP_D2_CN {
        sum += term.n * (term.t as f64) * delta.powi(term.d) * tau.powf((term.t - 1) as f64)
            * (-delta.powi(term.c)).exp();
    }

    // Gaussian terms (i=1 to i=3): ∂(Σᵢ nᵢδᵈⁱτᵗⁱexp[-αᵢ(δ-εᵢ)²-βᵢ(τ-γᵢ)²])/∂τ
    for term in &RES_GAUSS {
        let exp_term = gauss_exp(term, delta, tau);
        let deriv = (term.t as f64) / tau - 2.0 * (term.b as f64) * (tau - term.g);
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * exp_term * deriv;
    }

    // Non-analytic terms (i=1 to i=2): ∂(Σᵢ nᵢΔᵇⁱδF)/∂τ
    for term in &RES_NON_ANAL {
        let d_1 = delta - 1.0;
        let d_1_2 = d_1 * d_1;
        let tita = (1.0 - tau) + term.A * d_1_2.powf(0.5 / term.bt);
        let f_val = (-term.C as f64 * d_1_2 - term.D as f64 * (tau - 1.0).powi(2)).exp();
        let f_t = -2.0 * term.D as f64 * f_val * (tau - 1.0);

        let tita2 = tita * tita;
        let delta_val = tita2 + term.B * d_1_2.powf(term.a);

        let (delta_b, delta_bt) = if delta_val == 0.0 {
            (0.0, 0.0)
        } else {
            let db = delta_val.powf(term.b);
            let dbt = term.b * delta_val.powf(term.b - 1.0) * (-2.0 * tita);
            (db, dbt)
        };

        sum += term.n * delta * (delta_bt * f_val + delta_b * f_t);
    }

    sum
}

/// Compute second derivative ∂²φʳ/∂τ² for residual part
pub fn d2phi_residual_dtau2(delta: f64, tau: f64) -> f64 {
    let mut sum = 0.0f64;

    // Polynomial terms (i=1 to i=7): ∂²(Σᵢ nᵢ·δᵈⁱ·τᵗⁱ)/∂τ² = nᵢ·tᵢ·(tᵢ-1)·δᵈⁱ·τᵗⁱ⁻²
    for term in &RES_POLY_D1 {
        let t = term.t;
        sum += term.n * t * (t - 1.0) * delta.powi(term.d) * tau.powf(t - 2.0);
    }

    // Exponential terms (i=8 to i=22): c=1 ∂²(Σᵢ nᵢδᵈⁱτᵗⁱexp(-δ))/∂τ²
    for term in &RES_EXP_D2_C1 {
        let t = term.t as f64;
        sum += term.n * t * (t - 1.0) * delta.powi(term.d) * tau.powf(t - 2.0) * (-delta).exp();
    }

    // Exponential terms (i=23 to i=42): c=2 ∂²(Σᵢ nᵢδᵈⁱτᵗⁱexp(-δ²))/∂τ²
    let delta_2 = delta * delta;
    for term in &RES_EXP_D2_C2 {
        let t = term.t as f64;
        sum += term.n * t * (t - 1.0) * delta.powi(term.d) * tau.powf(t - 2.0) * (-delta_2).exp();
    }

    // Exponential terms (i=43 to i=51): ∂²(Σᵢ nᵢδᵈⁱτᵗⁱexp(-δᶜ))/∂τ²
    for term in &RES_EXP_D2_CN {
        let t = term.t as f64;
        sum += term.n * t * (t - 1.0) * delta.powi(term.d) * tau.powf(t - 2.0)
            * (-delta.powi(term.c)).exp();
    }

    // Gaussian terms (i=1 to i=3): ∂²(Σᵢ nᵢδᵈⁱτᵗⁱexp[-αᵢ(δ-εᵢ)²-βᵢ(τ-γᵢ)²])/∂τ²
    for term in &RES_GAUSS {
        let t_g = tau - term.g;
        let exp_term = gauss_exp(term, delta, tau);
        let t = term.t as f64;
        let b = term.b as f64;

        // Second derivative of τᵗ * exp[-β(τ-γ)²]
        // Let g(τ) = τᵗ * exp[-β(τ-γ)²]
        // g'(τ) = τᵗ * exp[-β(τ-γ)²] * (t/τ - 2β(τ-γ))
        // g''(τ) = τᵗ * exp[-β(τ-γ)²] * [(t/τ - 2β(τ-γ))² + (-t/τ² - 2β)]
        let factor = t / tau - 2.0 * b * t_g;
        let d2_factor = -t / (tau * tau) - 2.0 * b;
        let d2_tau_exp = tau.powf(t) * exp_term * (factor * factor + d2_factor);

        sum += term.n * delta.powi(term.d) * d2_tau_exp;
    }

    // Non-analytic terms: second derivative with respect to tau
    // φʳ_nonanal = nᵢ * Δᵇⁱ * δ * F(δ,τ)
    // where Δ = θ² + B*(δ-1)²ᵃ, θ = (1-τ) + A*(δ-1)²ᵇᵗ, F = exp[-C*(δ-1)² - D*(τ-1)²]
    for term in &RES_NON_ANAL {
        let d_1 = delta - 1.0;
        let d_1_2 = d_1 * d_1;
        let tita = (1.0 - tau) + term.A * d_1_2.powf(0.5 / term.bt);
        let f_val = (-term.C as f64 * d_1_2 - term.D as f64 * (tau - 1.0).powi(2)).exp();
        
        // First and second derivatives of F with respect to τ
        let f_t = -2.0 * term.D as f64 * f_val * (tau - 1.0);
        let f_tt = -2.0 * term.D as f64 * (f_val - 2.0 * term.D as f64 * f_val * (tau - 1.0) * (tau - 1.0));
    
        let tita2 = tita * tita;
        let delta_val = tita2 + term.B * d_1_2.powf(term.a);
        
        // First and second derivatives of Δ with respect to τ
        // ∂Δ/∂τ = 2*θ*∂θ/∂τ = 2*θ*(-1) = -2*θ
        let delta_t = -2.0 * tita;
        // ∂²Δ/∂τ² = 2*(∂θ/∂τ)² + 2*θ*∂²θ/∂τ² = 2*(-1)² + 0 = 2
        let delta_tt = 2.0;
    
        let (delta_b, delta_bt, delta_btt) = if delta_val == 0.0 {
            (0.0, 0.0, 0.0)
        } else {
            let db = delta_val.powf(term.b);
            let dbt = term.b * delta_val.powf(term.b - 1.0) * delta_t;
            let dbtt = term.b * (
                (term.b - 1.0) * delta_val.powf(term.b - 2.0) * delta_t * delta_t
                + delta_val.powf(term.b - 1.0) * delta_tt
            );
            (db, dbt, dbtt)
        };
    
        // Second derivative of the term nᵢ * (Δᵇⁱ * δ * F)
        // ∂²/∂τ²(Δᵇ * δ * F) = δ * [Δᵇᵗᵗ * F + 2*Δᵇᵗ * Fᵗ + Δᵇ * Fᵗᵗ]
        let term_tt = term.n * delta * (
            delta_btt * f_val
            + 2.0 * delta_bt * f_t
            + delta_b * f_tt
        );
        sum += term_tt;
    }

    sum
}

/// Compute mixed derivative ∂²φʳ/∂δ∂τ for residual part
pub fn d2phi_residual_ddelta_dtau(delta: f64, tau: f64) -> f64 {
    let mut sum = 0.0f64;

    // Polynomial terms (i=1 to i=7): ∂²(Σᵢ nᵢ·δᵈⁱ·τᵗⁱ)/∂δ∂τ = nᵢ·dᵢ·tᵢ·δᵈⁱ⁻¹·τᵗⁱ⁻¹
    for term in &RES_POLY_D1 {
        sum += term.n * (term.d as f64) * term.t * delta.powi(term.d - 1) * tau.powf(term.t - 1.0);
    }

    // Exponential terms (i=8 to i=22): c=1 ∂²(Σᵢ nᵢδᵈⁱτᵗⁱexp(-δ))/∂δ∂τ
    for term in &RES_EXP_D2_C1 {
        let deriv_delta = (term.d as f64) - delta;
        sum += term.n * (term.t as f64) * delta.powi(term.d - 1) * tau.powf((term.t - 1) as f64)
            * deriv_delta* (-delta).exp()
    }

    // Exponential terms (i=23 to i=42): c=2 
    let delta_2 = delta * delta;
    for term in &RES_EXP_D2_C2 {
        let deriv_delta = (term.d as f64) - 2.0* delta_2;
        sum += term.n * (term.t as f64) * delta.powi(term.d - 1) * tau.powf((term.t - 1) as f64)
            * deriv_delta * (-delta_2).exp();
    }

    // Exponential terms (i=43 to i=51): c=3,4,6
    for term in &RES_EXP_D2_CN {
        let delta_c = delta.powi(term.c);
        let deriv_delta = (term.d as f64) - (term.c as f64)* delta_c;
        sum += term.n * (term.t as f64) * delta.powi(term.d - 1) * tau.powf((term.t - 1) as f64)
            * deriv_delta * (-delta_c).exp();     
    }

    // Gaussian terms (i=1 to i=3): ∂²(Σᵢ nᵢδᵈⁱτᵗⁱexp[-αᵢ(δ-εᵢ)²-βᵢ(τ-γᵢ)²])/∂δ∂τ
    for term in &RES_GAUSS {
        let exp_term = gauss_exp(term, delta, tau);

        // Mixed derivative
        let deriv_delta = (term.d as f64) / delta - 2.0 * (term.a as f64) * (delta - term.e as f64);
        let deriv_tau = (term.t as f64) / tau - 2.0 * (term.b as f64) * (tau - term.g);

        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * exp_term * deriv_delta * deriv_tau;
    }

    // Non-analytic terms: mixed derivative
    // φʳ_nonanal = nᵢ * Δᵇⁱ * δ * F(δ,τ)
    // where Δ = θ² + B*(δ-1)²ᵃ, θ = (1-τ) + A*(δ-1)²ᵇᵗ, F = exp[-C*(δ-1)² - D*(τ-1)²]
    for term in &RES_NON_ANAL {
        let d_1 = delta - 1.0;
        let d_1_2 = d_1 * d_1;
        let tita = (1.0 - tau) + term.A * d_1_2.powf(0.5 / term.bt);
        let f_val = (-term.C as f64 * d_1_2 - term.D as f64 * (tau - 1.0).powi(2)).exp();
        
        // Derivatives of F
        let f_d = -2.0 * term.C as f64 * f_val * d_1;
        let f_t = -2.0 * term.D as f64 * f_val * (tau - 1.0);
        let f_dt = 4.0 * term.C as f64 * term.D as f64 * f_val * d_1 * (tau - 1.0);

        let tita2 = tita * tita;
        let delta_val = tita2 + term.B * d_1_2.powf(term.a);
        
        // Derivatives of Δ with respect to δ
        let tita_d = term.A / term.bt * d_1_2.powf(0.5 / term.bt - 1.0) * d_1;
        let delta_d = 2.0 * tita * tita_d + 2.0 * term.B * term.a * d_1_2.powf(term.a - 1.0) * d_1;
        
        // Derivatives of Δ with respect to τ
        let delta_t = -2.0 * tita;
        
        // Mixed derivative of Δ
        // ∂²Δ/∂δ∂τ = 2*(∂θ/∂δ)*(∂θ/∂τ) = 2*tita_d*(-1) = -2*tita_d
        let delta_dt = -2.0 * tita_d;

        let (delta_b, delta_bd, delta_bt, delta_bdt) = if delta_val == 0.0 {
            (0.0, 0.0, 0.0, 0.0)
        } else {
            let db = delta_val.powf(term.b);
            let dbd = term.b * delta_val.powf(term.b - 1.0) * delta_d;
            let dbt = term.b * delta_val.powf(term.b - 1.0) * delta_t;
            let dbdt = term.b * (
                (term.b - 1.0) * delta_val.powf(term.b - 2.0) * delta_d * delta_t
                + delta_val.powf(term.b - 1.0) * delta_dt
            );
            (db, dbd, dbt, dbdt)
        };

        // Mixed derivative of the term nᵢ * (Δᵇⁱ * δ * F)
        // G = Δᵇ * δ * F
        // ∂G/∂δ = Δᵇᵈ * δ * F + Δᵇ * F + Δᵇ * δ * F_d
        // ∂²G/∂δ∂τ = Δᵇᵈᵗ * δ * F + Δᵇᵈ * δ * F_t + Δᵇᵗ * F + Δᵇ * F_t + Δᵇᵗ * δ * F_d + Δᵇ * δ * F_dt
        let term_dt = term.n * (
            delta_bdt * delta * f_val
            + delta_bd * delta * f_t
            + delta_bt * f_val
            + delta_b * f_t
            + delta_bt * delta * f_d
            + delta_b * delta * f_dt
        );
        sum += term_dt;
    }

    sum
}
