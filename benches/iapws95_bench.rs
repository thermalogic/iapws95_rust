use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId, Throughput};
use iapws95::{
    tr2p, tr2u, tr2h, tr2s, tr2cv, tr2cp, tr2w, tr2jt, tr2itt, tr2beta_s,
    tr2p_direct, tr2u_direct, tr2h_direct, tr2s_direct, tr2cv_direct,
    tr2cp_direct, tr2w_direct, tr2jt_direct, tr2itt_direct, tr2beta_s_direct,
};

fn benchmark_single_properties(c: &mut Criterion) {
    let t_c = 25.0;
    let rho = 997.0;

    c.bench_function("tr2p (pressure)", |b| b.iter(|| tr2p(t_c, rho)));
    c.bench_function("tr2u (internal energy)", |b| b.iter(|| tr2u(t_c, rho)));
    c.bench_function("tr2h (enthalpy)", |b| b.iter(|| tr2h(t_c, rho)));
    c.bench_function("tr2s (entropy)", |b| b.iter(|| tr2s(t_c, rho)));
    c.bench_function("tr2cv (cv)", |b| b.iter(|| tr2cv(t_c, rho)));
    c.bench_function("tr2cp (cp)", |b| b.iter(|| tr2cp(t_c, rho)));
    c.bench_function("tr2w (speed of sound)", |b| b.iter(|| tr2w(t_c, rho)));
    c.bench_function("tr2jt (Joule-Thomson)", |b| b.iter(|| tr2jt(t_c, rho)));
    c.bench_function("tr2itt (isothermal throttling)", |b| b.iter(|| tr2itt(t_c, rho)));
    c.bench_function("tr2beta_s (isentropic temp-pressure)", |b| b.iter(|| tr2beta_s(t_c, rho)));
}

fn benchmark_steam_properties(c: &mut Criterion) {
    let t_c = 500.0;
    let rho = 5.0;

    c.bench_function("tr2p (steam)", |b| b.iter(|| tr2p(t_c, rho)));
    c.bench_function("tr2h (steam)", |b| b.iter(|| tr2h(t_c, rho)));
    c.bench_function("tr2cp (steam)", |b| b.iter(|| tr2cp(t_c, rho)));
    c.bench_function("tr2w (steam)", |b| b.iter(|| tr2w(t_c, rho)));
}

fn benchmark_multiple_states(c: &mut Criterion) {
    let states = [
        (25.0, 997.0),
        (100.0, 0.6),
        (200.0, 10.0),
        (350.0, 950.0),
        (500.0, 5.0),
    ];

    c.bench_function("all properties x 5 states", |b| {
        b.iter(|| {
            for (t_c, rho) in states {
                let _ = tr2p(t_c, rho);
                let _ = tr2u(t_c, rho);
                let _ = tr2h(t_c, rho);
                let _ = tr2s(t_c, rho);
                let _ = tr2cv(t_c, rho);
                let _ = tr2cp(t_c, rho);
                let _ = tr2w(t_c, rho);
                let _ = tr2jt(t_c, rho);
                let _ = tr2itt(t_c, rho);
                let _ = tr2beta_s(t_c, rho);
            }
        })
    });
}

// ========================================================================
// Precomputed vs Direct powi Performance Comparison
// ========================================================================

fn benchmark_precomputed_vs_direct(c: &mut Criterion) {
    let test_cases = [
        ("liquid", 25.0, 997.0),
        ("steam", 500.0, 5.0),
        ("supercritical", 400.0, 300.0),
        ("near_critical", 647.0, 322.0),
        ("high_temp", 1000.0, 10.0),
    ];

    let mut group = c.benchmark_group("precomputed_vs_direct");
    group.throughput(Throughput::Elements(1));

    for (name, t_c, rho) in test_cases {
        group.bench_with_input(BenchmarkId::new("pressure_precomputed", name), &(t_c, rho), |b, &(t_c, rho)| {
            b.iter(|| tr2p(t_c, rho))
        });
        group.bench_with_input(BenchmarkId::new("pressure_direct", name), &(t_c, rho), |b, &(t_c, rho)| {
            b.iter(|| tr2p_direct(t_c, rho))
        });

        group.bench_with_input(BenchmarkId::new("enthalpy_precomputed", name), &(t_c, rho), |b, &(t_c, rho)| {
            b.iter(|| tr2h(t_c, rho))
        });
        group.bench_with_input(BenchmarkId::new("enthalpy_direct", name), &(t_c, rho), |b, &(t_c, rho)| {
            b.iter(|| tr2h_direct(t_c, rho))
        });

        group.bench_with_input(BenchmarkId::new("entropy_precomputed", name), &(t_c, rho), |b, &(t_c, rho)| {
            b.iter(|| tr2s(t_c, rho))
        });
        group.bench_with_input(BenchmarkId::new("entropy_direct", name), &(t_c, rho), |b, &(t_c, rho)| {
            b.iter(|| tr2s_direct(t_c, rho))
        });

        group.bench_with_input(BenchmarkId::new("cp_precomputed", name), &(t_c, rho), |b, &(t_c, rho)| {
            b.iter(|| tr2cp(t_c, rho))
        });
        group.bench_with_input(BenchmarkId::new("cp_direct", name), &(t_c, rho), |b, &(t_c, rho)| {
            b.iter(|| tr2cp_direct(t_c, rho))
        });

        group.bench_with_input(BenchmarkId::new("speed_of_sound_precomputed", name), &(t_c, rho), |b, &(t_c, rho)| {
            b.iter(|| tr2w(t_c, rho))
        });
        group.bench_with_input(BenchmarkId::new("speed_of_sound_direct", name), &(t_c, rho), |b, &(t_c, rho)| {
            b.iter(|| tr2w_direct(t_c, rho))
        });
    }

    group.finish();
}

fn benchmark_all_properties_comparison(c: &mut Criterion) {
    let t_c = 25.0;
    let rho = 997.0;

    let mut group = c.benchmark_group("all_properties_comparison_liquid");

    group.bench_function("precomputed", |b| {
        b.iter(|| {
            let _ = tr2p(t_c, rho);
            let _ = tr2u(t_c, rho);
            let _ = tr2h(t_c, rho);
            let _ = tr2s(t_c, rho);
            let _ = tr2cv(t_c, rho);
            let _ = tr2cp(t_c, rho);
            let _ = tr2w(t_c, rho);
            let _ = tr2jt(t_c, rho);
            let _ = tr2itt(t_c, rho);
            let _ = tr2beta_s(t_c, rho);
        })
    });

    group.bench_function("direct", |b| {
        b.iter(|| {
            let _ = tr2p_direct(t_c, rho);
            let _ = tr2u_direct(t_c, rho);
            let _ = tr2h_direct(t_c, rho);
            let _ = tr2s_direct(t_c, rho);
            let _ = tr2cv_direct(t_c, rho);
            let _ = tr2cp_direct(t_c, rho);
            let _ = tr2w_direct(t_c, rho);
            let _ = tr2jt_direct(t_c, rho);
            let _ = tr2itt_direct(t_c, rho);
            let _ = tr2beta_s_direct(t_c, rho);
        })
    });

    group.finish();

    let t_c = 500.0;
    let rho = 5.0;

    let mut group = c.benchmark_group("all_properties_comparison_steam");

    group.bench_function("precomputed", |b| {
        b.iter(|| {
            let _ = tr2p(t_c, rho);
            let _ = tr2u(t_c, rho);
            let _ = tr2h(t_c, rho);
            let _ = tr2s(t_c, rho);
            let _ = tr2cv(t_c, rho);
            let _ = tr2cp(t_c, rho);
            let _ = tr2w(t_c, rho);
            let _ = tr2jt(t_c, rho);
            let _ = tr2itt(t_c, rho);
            let _ = tr2beta_s(t_c, rho);
        })
    });

    group.bench_function("direct", |b| {
        b.iter(|| {
            let _ = tr2p_direct(t_c, rho);
            let _ = tr2u_direct(t_c, rho);
            let _ = tr2h_direct(t_c, rho);
            let _ = tr2s_direct(t_c, rho);
            let _ = tr2cv_direct(t_c, rho);
            let _ = tr2cp_direct(t_c, rho);
            let _ = tr2w_direct(t_c, rho);
            let _ = tr2jt_direct(t_c, rho);
            let _ = tr2itt_direct(t_c, rho);
            let _ = tr2beta_s_direct(t_c, rho);
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    benchmark_single_properties,
    benchmark_steam_properties,
    benchmark_multiple_states,
    benchmark_precomputed_vs_direct,
    benchmark_all_properties_comparison,
);
criterion_main!(benches);
