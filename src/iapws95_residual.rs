//! IAPWS-95 Residual Part Implementation (Internal Module)
//!
//! **This module is internal (`pub(crate)`) and not exposed to external users.**
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
    a: i32,
    b: i32,
    g: f64,
    e: i32,
}

/// Non-analytic term: (alpha, beta, B, n, C, D, A, bt)
struct NonAnalTerm {
    a: f64,
    b: f64,
    B: f64,
    n: f64,
    C: i32,
    D: i32,
    A: f64,
    bt: f64,
}

const RES_POLY_D1: [PolyTerm; 7] = [
    PolyTerm { d: 1, t: -0.5, n: 0.12533547935523e-1 },
    PolyTerm { d: 1, t:  0.875, n: 0.78957634722828e1 },
    PolyTerm { d: 1, t:  1.0, n: -0.87803203303561e1 },
    PolyTerm { d: 2, t:  0.5, n: 0.31802509345418 },
    PolyTerm { d: 2, t:  0.75, n: -0.26145533859358 },
    PolyTerm { d: 3, t:  0.375, n: -0.78199751687981e-2 },
    PolyTerm { d: 4, t:  1.0, n: 0.88089493102134e-2 },
];

const RES_EXP_D2_C1: [ExpTermC1C2; 15] = [
    ExpTermC1C2 { d: 1, t:   4, n: -0.66856572307965 },
    ExpTermC1C2 { d: 1, t:   6, n:  0.20433810950965 },
    ExpTermC1C2 { d: 1, t:  12, n: -0.66212605039687e-4 },
    ExpTermC1C2 { d: 2, t:   1, n: -0.19232721156002 },
    ExpTermC1C2 { d: 2, t:   5, n: -0.25709043003438 },
    ExpTermC1C2 { d: 3, t:   4, n:  0.16074868486251 },
    ExpTermC1C2 { d: 4, t:   2, n: -0.40092828925807e-1 },
    ExpTermC1C2 { d: 4, t:  13, n:  0.39343422603254e-6 },
    ExpTermC1C2 { d: 5, t:   9, n: -0.75941377088144e-5 },
    ExpTermC1C2 { d: 7, t:   3, n:  0.56250979351888e-3 },
    ExpTermC1C2 { d: 9, t:   4, n: -0.15608652257135e-4 },
    ExpTermC1C2 { d: 10, t: 11, n:  0.11537996422951e-8 },
    ExpTermC1C2 { d: 11, t:  4, n:  0.36582165144204e-6 },
    ExpTermC1C2 { d: 13, t: 13, n: -0.13251180074668e-11 },
    ExpTermC1C2 { d: 15, t:  1, n: -0.62639586912454e-9 },
];

const RES_EXP_D2_C2: [ExpTermC1C2; 20] = [
    ExpTermC1C2 { d: 1, t:   7, n: -0.10793600908932 },
    ExpTermC1C2 { d: 2, t:   1, n:  0.17611491008752e-1 },
    ExpTermC1C2 { d: 2, t:   9, n:  0.22132295167546 },
    ExpTermC1C2 { d: 2, t:  10, n: -0.40247669763528 },
    ExpTermC1C2 { d: 3, t:  10, n:  0.58083399985759 },
    ExpTermC1C2 { d: 4, t:   3, n:  0.49969146990806e-2 },
    ExpTermC1C2 { d: 4, t:   7, n: -0.31358700712549e-1 },
    ExpTermC1C2 { d: 4, t:  10, n: -0.74315929710341 },
    ExpTermC1C2 { d: 5, t:  10, n:  0.47807329915480 },
    ExpTermC1C2 { d: 6, t:   6, n:  0.20527940895948e-1 },
    ExpTermC1C2 { d: 6, t:  10, n: -0.13636435110343 },
    ExpTermC1C2 { d: 7, t:  10, n:  0.14180634400617e-1 },
    ExpTermC1C2 { d: 9, t:   1, n:  0.83326504880713e-2 },
    ExpTermC1C2 { d: 9, t:   2, n: -0.29052336009585e-1 },
    ExpTermC1C2 { d: 9, t:   3, n:  0.38615085574206e-1 },
    ExpTermC1C2 { d: 9, t:   4, n: -0.20393486513704e-1 },
    ExpTermC1C2 { d: 9, t:   8, n: -0.16554050063734e-2 },
    ExpTermC1C2 { d: 10, t:  6, n:  0.19955571979541e-2 },
    ExpTermC1C2 { d: 10, t:  9, n:  0.15870308324157e-3 },
    ExpTermC1C2 { d: 12, t:  8, n: -0.16388568342530e-4 },
];

const RES_EXP_D2_CN: [ExpTermCN; 9] = [
    ExpTermCN { c: 3, d: 3, t: 16, n:  0.43613615723811e-1 },
    ExpTermCN { c: 3, d: 4, t: 22, n:  0.34994005463765e-1 },
    ExpTermCN { c: 3, d: 4, t: 23, n: -0.76788197844621e-1 },
    ExpTermCN { c: 3, d: 5, t: 23, n:  0.22446277332006e-1 },
    ExpTermCN { c: 4, d: 14, t: 10, n: -0.62689710414685e-4 },
    ExpTermCN { c: 6, d: 3, t: 50, n: -0.55711118565645e-9 },
    ExpTermCN { c: 6, d: 6, t: 44, n: -0.19905718354408 },
    ExpTermCN { c: 6, d: 6, t: 46, n:  0.31777497330738 },
    ExpTermCN { c: 6, d: 6, t: 50, n: -0.11841182425981 },
];

const RES_GAUSS: [GaussTerm; 3] = [
    GaussTerm { d: 3, t: 0, n: -0.31306260323435e2, a: 20, b: 150, g: 1.21, e: 1 },
    GaussTerm { d: 3, t: 1, n:  0.31546140237781e2, a: 20, b: 150, g: 1.21, e: 1 },
    GaussTerm { d: 3, t: 4, n: -0.25213154341695e4, a: 20, b: 250, g: 1.25, e: 1 },
];

const RES_NON_ANAL: [NonAnalTerm; 2] = [
    NonAnalTerm { a: 3.5, b: 0.85, B: 0.2, n: -0.14874640856724, C: 28, D: 700, A: 0.32, bt: 0.3 },
    NonAnalTerm { a: 3.5, b: 0.95, B: 0.2, n:  0.31806110878444, C: 32, D: 800, A: 0.32, bt: 0.3 },
];

#[inline]
fn gauss_exp(term: &GaussTerm, delta: f64, tau: f64) -> f64 {
    let d_e = delta - term.e as f64;
    let t_g = tau - term.g;
    (-(term.a as f64) * d_e * d_e - (term.b as f64) * t_g * t_g).exp()
}

#[inline]
pub fn phi_residual(delta: f64, tau: f64) -> f64 {
    let mut sum = 0.0f64;

    let delta_2 = delta * delta;
    let delta_3 = delta_2 * delta;
    let delta_4 = delta_3 * delta;
    let delta_6 = delta_4 * delta_2;

    for term in &RES_POLY_D1 {
        sum += term.n * delta.powi(term.d) * tau.powf(term.t);
    }

    let exp_delta = (-delta).exp();
    for term in &RES_EXP_D2_C1 {
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * exp_delta;
    }

    let exp_delta_2 = (-delta_2).exp();
    for term in &RES_EXP_D2_C2 {
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * exp_delta_2;
    }

    let exp_delta_3 = (-delta_3).exp();
    let exp_delta_4 = (-delta_4).exp();
    let exp_delta_6 = (-delta_6).exp();

    for term in &RES_EXP_D2_CN[0..4] {
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * exp_delta_3;
    }
    let term = &RES_EXP_D2_CN[4];
    sum += term.n * delta.powi(term.d) * tau.powi(term.t) * exp_delta_4;
    for term in &RES_EXP_D2_CN[5..9] {
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * exp_delta_6;
    }

    for term in &RES_GAUSS {
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * gauss_exp(term, delta, tau);
    }

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

#[inline]
pub fn dphi_residual_ddelta(delta: f64, tau: f64) -> f64 {
    let mut sum = 0.0f64;

    let delta_2 = delta * delta;
    let delta_3 = delta_2 * delta;
    let delta_4 = delta_3 * delta;
    let delta_6 = delta_4 * delta_2;

    for term in &RES_POLY_D1 {
        sum += term.n * (term.d as f64) * delta.powi(term.d - 1) * tau.powf(term.t);
    }

    let exp_delta = (-delta).exp();
    for term in &RES_EXP_D2_C1 {
        let deriv = (term.d as f64) - delta;
        sum += term.n * exp_delta * delta.powi(term.d - 1) * tau.powi(term.t) * deriv;
    }

    let exp_delta_2 = (-delta_2).exp();
    for term in &RES_EXP_D2_C2 {
        let deriv = (term.d as f64) - 2.0 * delta_2;
        sum += term.n * exp_delta_2 * delta.powi(term.d - 1) * tau.powi(term.t) * deriv;
    }

    let exp_delta_3 = (-delta_3).exp();
    let exp_delta_4 = (-delta_4).exp();
    let exp_delta_6 = (-delta_6).exp();
    for term in &RES_EXP_D2_CN[0..4] {
        let deriv = (term.d as f64) - 3.0 * delta_3;
        sum += term.n * exp_delta_3 * delta.powi(term.d - 1) * tau.powi(term.t) * deriv;
    }
    let term = &RES_EXP_D2_CN[4];
    let deriv = (term.d as f64) - 4.0 * delta_4;
    sum += term.n * exp_delta_4 * delta.powi(term.d - 1) * tau.powi(term.t) * deriv;
    for term in &RES_EXP_D2_CN[5..9] {
        let deriv = (term.d as f64) - 6.0 * delta_6;
        sum += term.n * exp_delta_6 * delta.powi(term.d - 1) * tau.powi(term.t) * deriv;
    }

    for term in &RES_GAUSS {
        let exp_term = gauss_exp(term, delta, tau);
        let deriv = (term.d as f64) / delta - 2.0 * (term.a as f64) * (delta - term.e as f64);
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * exp_term * deriv;
    }

    for term in &RES_NON_ANAL {
        let d_1 = delta - 1.0;
        let d_1_2 = d_1 * d_1;
        let tita = (1.0 - tau) + term.A * d_1_2.powf(0.5 / term.bt);
        let f_val = (-term.C as f64 * d_1_2 - term.D as f64 * (tau - 1.0) * (tau - 1.0)).exp();
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

#[inline]
pub fn d2phi_residual_ddelta2(delta: f64, tau: f64) -> f64 {
    let mut sum = 0.0f64;

    let delta_2 = delta * delta;
    let delta_3 = delta_2 * delta;
    let delta_4 = delta_3 * delta;
    let delta_6 = delta_4 * delta_2;

    for term in &RES_POLY_D1 {
        let d = term.d as f64;
        if term.d >= 2 {
            sum += term.n * d * (d - 1.0) * delta.powi(term.d - 2) * tau.powf(term.t);
        }
    }

    let exp_delta = (-delta).exp();
    for term in &RES_EXP_D2_C1 {
        let d = term.d as f64;
        let delta_d_minus_2 = delta.powi(term.d - 2);
        let d_term = delta_d_minus_2 * tau.powi(term.t) * ((d - delta) * (d - 1.0 - delta) - delta);
        sum += term.n * exp_delta * d_term;
    }

    let exp_delta_2 = (-delta_2).exp();
    for term in &RES_EXP_D2_C2 {
        let d = term.d as f64;
        let delta_d_minus_2 = delta.powi(term.d - 2);
        let d_term = delta_d_minus_2 * tau.powi(term.t) * ((d - 1.0 - 2.0 * delta_2) * (d - 2.0 * delta_2) - 4.0 * delta_2);
        sum += term.n * exp_delta_2 * d_term;
    }

    let exp_delta_3 = (-delta_3).exp();
    let exp_delta_4 = (-delta_4).exp();
    let exp_delta_6 = (-delta_6).exp();
    for term in &RES_EXP_D2_CN[0..4] {
        let d = term.d as f64;
        let delta_d_minus_2 = delta.powi(term.d - 2);
        let d_term = delta_d_minus_2 * tau.powi(term.t) * ((d - 1.0 - 3.0 * delta_3) * (d - 3.0 * delta_3) - 9.0 * delta_3);
        sum += term.n * exp_delta_3 * d_term;
    }
    let term = &RES_EXP_D2_CN[4];
    let d = term.d as f64;
    let delta_d_minus_2 = delta.powi(term.d - 2);
    let d_term = delta_d_minus_2 * tau.powi(term.t) * ((d - 1.0 - 4.0 * delta_4) * (d - 4.0 * delta_4) - 16.0 * delta_4);
    sum += term.n * exp_delta_4 * d_term;
    for term in &RES_EXP_D2_CN[5..9] {
        let d = term.d as f64;
        let delta_d_minus_2 = delta.powi(term.d - 2);
        let d_term = delta_d_minus_2 * tau.powi(term.t) * ((d - 1.0 - 6.0 * delta_6) * (d - 6.0 * delta_6) - 36.0 * delta_6);
        sum += term.n * exp_delta_6 * d_term;
    }

    for term in &RES_GAUSS {
        let d = term.d as f64;
        let a = term.a as f64;
        let e = term.e as f64;
        let delta_e = delta - e;
        let exp_term = gauss_exp(term, delta, tau);

        let first_deriv_factor = d / delta - 2.0 * a * delta_e;
        let second_deriv_correction = -d / delta_2 - 2.0 * a;
        let factor = first_deriv_factor * first_deriv_factor + second_deriv_correction;

        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * exp_term * factor;
    }

    for term in &RES_NON_ANAL {
        let d_1 = delta - 1.0;
        let d_1_2 = d_1 * d_1;
        let tita = (1.0 - tau) + term.A * d_1_2.powf(0.5 / term.bt);
        let f_val = (-term.C as f64 * d_1_2 - term.D as f64 * (tau - 1.0) * (tau - 1.0)).exp();

        let f_d = -2.0 * term.C as f64 * f_val * d_1;
        let f_dd = -2.0 * term.C as f64 * (f_val + d_1 * f_d);

        let tita2 = tita * tita;
        let delta_val = tita2 + term.B * d_1_2.powf(term.a);

        let delta_d = d_1 * (
            term.A * tita * 2.0 / term.bt * d_1_2.powf(0.5 / term.bt - 1.0)
            + 2.0 * term.B * term.a * d_1_2.powf(term.a - 1.0)
        );

        let tita_d = term.A * d_1 * 1.0 / term.bt * d_1_2.powf(0.5 / term.bt - 1.0);
        let tita_dd = term.A * (1.0 / term.bt) * (
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

        let term_dd = term.n * (
            delta_bdd * delta * f_val
            + 2.0 * delta_bd * (f_val + delta * f_d)
            + delta_b * (2.0 * f_d + delta * f_dd)
        );
        sum += term_dd;
    }

    sum
}

#[inline]
pub fn dphi_residual_dtau(delta: f64, tau: f64) -> f64 {
    let mut sum = 0.0f64;

    let delta_2 = delta * delta;
    let delta_3 = delta_2 * delta;
    let delta_4 = delta_3 * delta;
    let delta_6 = delta_4 * delta_2;

    for term in &RES_POLY_D1 {
        sum += term.n * term.t * delta.powi(term.d) * tau.powf(term.t - 1.0);
    }

    let exp_delta = (-delta).exp();
    for term in &RES_EXP_D2_C1 {
        sum += term.n * (term.t as f64) * delta.powi(term.d) * tau.powi(term.t - 1) * exp_delta;
    }

    let exp_delta_2 = (-delta_2).exp();
    for term in &RES_EXP_D2_C2 {
        sum += term.n * (term.t as f64) * delta.powi(term.d) * tau.powi(term.t - 1) * exp_delta_2;
    }

    let exp_delta_3 = (-delta_3).exp();
    let exp_delta_4 = (-delta_4).exp();
    let exp_delta_6 = (-delta_6).exp();
    for term in &RES_EXP_D2_CN[0..4] {
        sum += term.n * (term.t as f64) * delta.powi(term.d) * tau.powi(term.t - 1) * exp_delta_3;
    }
    let term = &RES_EXP_D2_CN[4];
    sum += term.n * (term.t as f64) * delta.powi(term.d) * tau.powi(term.t - 1) * exp_delta_4;
    for term in &RES_EXP_D2_CN[5..9] {
        sum += term.n * (term.t as f64) * delta.powi(term.d) * tau.powi(term.t - 1) * exp_delta_6;
    }

    for term in &RES_GAUSS {
        let exp_term = gauss_exp(term, delta, tau);
        let deriv = (term.t as f64) / tau - 2.0 * (term.b as f64) * (tau - term.g);
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * exp_term * deriv;
    }

    for term in &RES_NON_ANAL {
        let d_1 = delta - 1.0;
        let d_1_2 = d_1 * d_1;
        let tita = (1.0 - tau) + term.A * d_1_2.powf(0.5 / term.bt);
        let f_val = (-term.C as f64 * d_1_2 - term.D as f64 * (tau - 1.0) * (tau - 1.0)).exp();
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

#[inline]
pub fn d2phi_residual_dtau2(delta: f64, tau: f64) -> f64 {
    let mut sum = 0.0f64;

    let delta_2 = delta * delta;
    let delta_3 = delta_2 * delta;
    let delta_4 = delta_3 * delta;
    let delta_6 = delta_4 * delta_2;

    for term in &RES_POLY_D1 {
        let t = term.t;
        sum += term.n * t * (t - 1.0) * delta.powi(term.d) * tau.powf(t - 2.0);
    }

    let exp_delta = (-delta).exp();
    for term in &RES_EXP_D2_C1 {
        let t = term.t as f64;
        if term.t >= 2 {
            sum += term.n * t * (t - 1.0) * delta.powi(term.d) * tau.powi(term.t - 2) * exp_delta;
        }
    }

    let exp_delta_2 = (-delta_2).exp();
    for term in &RES_EXP_D2_C2 {
        let t = term.t as f64;
        if term.t >= 2 {
            sum += term.n * t * (t - 1.0) * delta.powi(term.d) * tau.powi(term.t - 2) * exp_delta_2;
        }
    }

    let exp_delta_3 = (-delta_3).exp();
    let exp_delta_4 = (-delta_4).exp();
    let exp_delta_6 = (-delta_6).exp();
    for term in &RES_EXP_D2_CN[0..4] {
        let t = term.t as f64;
        if term.t >= 2 {
            sum += term.n * t * (t - 1.0) * delta.powi(term.d) * tau.powi(term.t - 2) * exp_delta_3;
        }
    }
    let term = &RES_EXP_D2_CN[4];
    let t = term.t as f64;
    if term.t >= 2 {
        sum += term.n * t * (t - 1.0) * delta.powi(term.d) * tau.powi(term.t - 2) * exp_delta_4;
    }
    for term in &RES_EXP_D2_CN[5..9] {
        let t = term.t as f64;
        if term.t >= 2 {
            sum += term.n * t * (t - 1.0) * delta.powi(term.d) * tau.powi(term.t - 2) * exp_delta_6;
        }
    }

    for term in &RES_GAUSS {
        let t_g = tau - term.g;
        let exp_term = gauss_exp(term, delta, tau);
        let t = term.t as f64;
        let b = term.b as f64;

        let factor = t / tau - 2.0 * b * t_g;
        let tau_2 = tau * tau;
        let d2_factor = -t / tau_2 - 2.0 * b;
        let d2_tau_exp = tau.powi(term.t) * exp_term * (factor * factor + d2_factor);

        sum += term.n * delta.powi(term.d) * d2_tau_exp;
    }

    for term in &RES_NON_ANAL {
        let d_1 = delta - 1.0;
        let d_1_2 = d_1 * d_1;
        let tita = (1.0 - tau) + term.A * d_1_2.powf(0.5 / term.bt);
        let f_val = (-term.C as f64 * d_1_2 - term.D as f64 * (tau - 1.0) * (tau - 1.0)).exp();

        let f_t = -2.0 * term.D as f64 * f_val * (tau - 1.0);
        let f_tt = -2.0 * term.D as f64 * (f_val - 2.0 * term.D as f64 * f_val * (tau - 1.0) * (tau - 1.0));

        let tita2 = tita * tita;
        let delta_val = tita2 + term.B * d_1_2.powf(term.a);

        let delta_t = -2.0 * tita;
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

        let term_tt = term.n * delta * (
            delta_btt * f_val
            + 2.0 * delta_bt * f_t
            + delta_b * f_tt
        );
        sum += term_tt;
    }

    sum
}

#[inline]
pub fn d2phi_residual_ddelta_dtau(delta: f64, tau: f64) -> f64 {
    let mut sum = 0.0f64;

    let delta_2 = delta * delta;
    let delta_3 = delta_2 * delta;
    let delta_4 = delta_3 * delta;
    let delta_6 = delta_4 * delta_2;

    for term in &RES_POLY_D1 {
        sum += term.n * (term.d as f64) * term.t * delta.powi(term.d - 1) * tau.powf(term.t - 1.0);
    }

    let exp_delta = (-delta).exp();
    for term in &RES_EXP_D2_C1 {
        let deriv_delta = (term.d as f64) - delta;
        sum += term.n * (term.t as f64) * delta.powi(term.d - 1) * tau.powi(term.t - 1)
            * deriv_delta * exp_delta;
    }

    let exp_delta_2 = (-delta_2).exp();
    for term in &RES_EXP_D2_C2 {
        let deriv_delta = (term.d as f64) - 2.0 * delta_2;
        sum += term.n * (term.t as f64) * delta.powi(term.d - 1) * tau.powi(term.t - 1)
            * deriv_delta * exp_delta_2;
    }

    let exp_delta_3 = (-delta_3).exp();
    let exp_delta_4 = (-delta_4).exp();
    let exp_delta_6 = (-delta_6).exp();
    for term in &RES_EXP_D2_CN[0..4] {
        let deriv_delta = (term.d as f64) - 3.0 * delta_3;
        sum += term.n * (term.t as f64) * delta.powi(term.d - 1) * tau.powi(term.t - 1)
            * deriv_delta * exp_delta_3;
    }
    let term = &RES_EXP_D2_CN[4];
    let deriv_delta = (term.d as f64) - 4.0 * delta_4;
    sum += term.n * (term.t as f64) * delta.powi(term.d - 1) * tau.powi(term.t - 1)
        * deriv_delta * exp_delta_4;
    for term in &RES_EXP_D2_CN[5..9] {
        let deriv_delta = (term.d as f64) - 6.0 * delta_6;
        sum += term.n * (term.t as f64) * delta.powi(term.d - 1) * tau.powi(term.t - 1)
            * deriv_delta * exp_delta_6;
    }

    for term in &RES_GAUSS {
        let exp_term = gauss_exp(term, delta, tau);
        let deriv_delta = (term.d as f64) / delta - 2.0 * (term.a as f64) * (delta - term.e as f64);
        let deriv_tau = (term.t as f64) / tau - 2.0 * (term.b as f64) * (tau - term.g);
        sum += term.n * delta.powi(term.d) * tau.powi(term.t) * exp_term * deriv_delta * deriv_tau;
    }

    for term in &RES_NON_ANAL {
        let d_1 = delta - 1.0;
        let d_1_2 = d_1 * d_1;
        let tita = (1.0 - tau) + term.A * d_1_2.powf(0.5 / term.bt);
        let f_val = (-term.C as f64 * d_1_2 - term.D as f64 * (tau - 1.0) * (tau - 1.0)).exp();

        let f_d = -2.0 * term.C as f64 * f_val * d_1;
        let f_t = -2.0 * term.D as f64 * f_val * (tau - 1.0);
        let f_dt = 4.0 * term.C as f64 * term.D as f64 * f_val * d_1 * (tau - 1.0);

        let tita2 = tita * tita;
        let delta_val = tita2 + term.B * d_1_2.powf(term.a);

        let tita_d = term.A / term.bt * d_1_2.powf(0.5 / term.bt - 1.0) * d_1;
        let delta_d = 2.0 * tita * tita_d + 2.0 * term.B * term.a * d_1_2.powf(term.a - 1.0) * d_1;
        let delta_t = -2.0 * tita;
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
