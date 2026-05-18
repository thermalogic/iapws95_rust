use std::time::Instant;
use iapws95::{
    tr2p, tr2h, tr2s, tr2cp, tr2w,
    tr2p_direct, tr2h_direct, tr2s_direct, tr2cp_direct, tr2w_direct,
};

const ITERATIONS: usize = 100_000;

fn main() {
    println!("=== IAPWS-95 Precomputed vs Direct powi Performance Test ===\n");

    let test_cases = [
        ("Liquid water (25°C, 997 kg/m³)", 25.0, 997.0),
        ("Steam (500°C, 5 kg/m³)", 500.0, 5.0),
        ("Supercritical (400°C, 300 kg/m³)", 400.0, 300.0),
    ];

    for (name, t_c, rho) in test_cases {
        println!("--- {} ---", name);

        let mut sum_p = 0.0;
        let start = Instant::now();
        for _ in 0..ITERATIONS { sum_p += tr2p(t_c, rho); }
        let ms_p = start.elapsed().as_secs_f64() * 1000.0;

        let mut sum_pd = 0.0;
        let start = Instant::now();
        for _ in 0..ITERATIONS { sum_pd += tr2p_direct(t_c, rho); }
        let ms_pd = start.elapsed().as_secs_f64() * 1000.0;

        println!("  Pressure:       precomputed={:.2}ms, direct={:.2}ms, ratio={:.2}x", ms_p, ms_pd, ms_pd/ms_p);

        let mut sum_h = 0.0;
        let start = Instant::now();
        for _ in 0..ITERATIONS { sum_h += tr2h(t_c, rho); }
        let ms_h = start.elapsed().as_secs_f64() * 1000.0;

        let mut sum_hd = 0.0;
        let start = Instant::now();
        for _ in 0..ITERATIONS { sum_hd += tr2h_direct(t_c, rho); }
        let ms_hd = start.elapsed().as_secs_f64() * 1000.0;

        println!("  Enthalpy:       precomputed={:.2}ms, direct={:.2}ms, ratio={:.2}x", ms_h, ms_hd, ms_hd/ms_h);

        let mut sum_s = 0.0;
        let start = Instant::now();
        for _ in 0..ITERATIONS { sum_s += tr2s(t_c, rho); }
        let ms_s = start.elapsed().as_secs_f64() * 1000.0;

        let mut sum_sd = 0.0;
        let start = Instant::now();
        for _ in 0..ITERATIONS { sum_sd += tr2s_direct(t_c, rho); }
        let ms_sd = start.elapsed().as_secs_f64() * 1000.0;

        println!("  Entropy:        precomputed={:.2}ms, direct={:.2}ms, ratio={:.2}x", ms_s, ms_sd, ms_sd/ms_s);

        let mut sum_cp = 0.0;
        let start = Instant::now();
        for _ in 0..ITERATIONS { sum_cp += tr2cp(t_c, rho); }
        let ms_cp = start.elapsed().as_secs_f64() * 1000.0;

        let mut sum_cpd = 0.0;
        let start = Instant::now();
        for _ in 0..ITERATIONS { sum_cpd += tr2cp_direct(t_c, rho); }
        let ms_cpd = start.elapsed().as_secs_f64() * 1000.0;

        println!("  Cp:             precomputed={:.2}ms, direct={:.2}ms, ratio={:.2}x", ms_cp, ms_cpd, ms_cpd/ms_cp);

        let mut sum_w = 0.0;
        let start = Instant::now();
        for _ in 0..ITERATIONS { sum_w += tr2w(t_c, rho); }
        let ms_w = start.elapsed().as_secs_f64() * 1000.0;

        let mut sum_wd = 0.0;
        let start = Instant::now();
        for _ in 0..ITERATIONS { sum_wd += tr2w_direct(t_c, rho); }
        let ms_wd = start.elapsed().as_secs_f64() * 1000.0;

        println!("  Speed of sound: precomputed={:.2}ms, direct={:.2}ms, ratio={:.2}x", ms_w, ms_wd, ms_wd/ms_w);

        println!("  Value check: p={:.6} vs {:.6}, h={:.6} vs {:.6}", sum_p, sum_pd, sum_h, sum_hd);
        println!();
    }
}
