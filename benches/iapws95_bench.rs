use criterion::{criterion_group, criterion_main, Criterion};
use iapws95::{tr2p, tr2u, tr2h, tr2s, tr2cv, tr2cp, tr2w, tr2jt, tr2itt, tr2beta_s};

fn benchmark_properties(c: &mut Criterion) {
    let states = [
        ("liquid", 25.0, 997.0),
        ("steam", 500.0, 5.0),
        ("supercritical", 400.0, 300.0),
    ];

    let mut group = c.benchmark_group("properties");

    for (name, t_c, rho) in states {
        group.bench_function(format!("{name}/pressure"), |b| b.iter(|| tr2p(t_c, rho)));
        group.bench_function(format!("{name}/energy"), |b| b.iter(|| tr2u(t_c, rho)));
        group.bench_function(format!("{name}/enthalpy"), |b| b.iter(|| tr2h(t_c, rho)));
        group.bench_function(format!("{name}/entropy"), |b| b.iter(|| tr2s(t_c, rho)));
        group.bench_function(format!("{name}/cv"), |b| b.iter(|| tr2cv(t_c, rho)));
        group.bench_function(format!("{name}/cp"), |b| b.iter(|| tr2cp(t_c, rho)));
        group.bench_function(format!("{name}/speed_of_sound"), |b| b.iter(|| tr2w(t_c, rho)));
        group.bench_function(format!("{name}/joule_thomson"), |b| b.iter(|| tr2jt(t_c, rho)));
        group.bench_function(format!("{name}/isothermal_throttling"), |b| b.iter(|| tr2itt(t_c, rho)));
        group.bench_function(format!("{name}/isentropic_temp_pressure"), |b| b.iter(|| tr2beta_s(t_c, rho)));
    }

    group.finish();
}

fn benchmark_all_properties(c: &mut Criterion) {
    let states = [
        (25.0, 997.0),
        (100.0, 0.6),
        (200.0, 10.0),
        (350.0, 950.0),
        (500.0, 5.0),
    ];

    c.bench_function("all_properties_x_5_states", |b| {
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

criterion_group!(
    benches,
    benchmark_properties,
    benchmark_all_properties
);
criterion_main!(benches);
