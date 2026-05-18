use criterion::{black_box, criterion_group, criterion_main, Criterion};
use iapws95::{tr2p, tr2u, tr2h, tr2s, tr2cv, tr2cp, tr2w, tr2jt, tr2itt, tr2beta_s};

fn benchmark_single_properties(c: &mut Criterion) {
    let t_c = 25.0;
    let rho = 997.0;

    c.bench_function("tr2p (pressure)", |b| b.iter(|| tr2p(black_box(t_c), black_box(rho))));
    c.bench_function("tr2u (internal energy)", |b| b.iter(|| tr2u(black_box(t_c), black_box(rho))));
    c.bench_function("tr2h (enthalpy)", |b| b.iter(|| tr2h(black_box(t_c), black_box(rho))));
    c.bench_function("tr2s (entropy)", |b| b.iter(|| tr2s(black_box(t_c), black_box(rho))));
    c.bench_function("tr2cv (cv)", |b| b.iter(|| tr2cv(black_box(t_c), black_box(rho))));
    c.bench_function("tr2cp (cp)", |b| b.iter(|| tr2cp(black_box(t_c), black_box(rho))));
    c.bench_function("tr2w (speed of sound)", |b| b.iter(|| tr2w(black_box(t_c), black_box(rho))));
    c.bench_function("tr2jt (Joule-Thomson)", |b| b.iter(|| tr2jt(black_box(t_c), black_box(rho))));
    c.bench_function("tr2itt (isothermal throttling)", |b| b.iter(|| tr2itt(black_box(t_c), black_box(rho))));
    c.bench_function("tr2beta_s (isentropic temp-pressure)", |b| b.iter(|| tr2beta_s(black_box(t_c), black_box(rho))));
}

fn benchmark_steam_properties(c: &mut Criterion) {
    let t_c = 500.0;
    let rho = 5.0;

    c.bench_function("tr2p (steam)", |b| b.iter(|| tr2p(black_box(t_c), black_box(rho))));
    c.bench_function("tr2h (steam)", |b| b.iter(|| tr2h(black_box(t_c), black_box(rho))));
    c.bench_function("tr2cp (steam)", |b| b.iter(|| tr2cp(black_box(t_c), black_box(rho))));
    c.bench_function("tr2w (steam)", |b| b.iter(|| tr2w(black_box(t_c), black_box(rho))));
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
                let _ = tr2p(black_box(t_c), black_box(rho));
                let _ = tr2u(black_box(t_c), black_box(rho));
                let _ = tr2h(black_box(t_c), black_box(rho));
                let _ = tr2s(black_box(t_c), black_box(rho));
                let _ = tr2cv(black_box(t_c), black_box(rho));
                let _ = tr2cp(black_box(t_c), black_box(rho));
                let _ = tr2w(black_box(t_c), black_box(rho));
                let _ = tr2jt(black_box(t_c), black_box(rho));
                let _ = tr2itt(black_box(t_c), black_box(rho));
                let _ = tr2beta_s(black_box(t_c), black_box(rho));
            }
        })
    });
}

criterion_group!(
    benches,
    benchmark_single_properties,
    benchmark_steam_properties,
    benchmark_multiple_states
);
criterion_main!(benches);
