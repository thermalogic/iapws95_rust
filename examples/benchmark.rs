use std::hint::black_box;
use std::time::Instant;
use iapws95::iapws95::{tr2p, tr2u, tr2h, tr2s, tr2cv, tr2cp, tr2w, tr2jt, tr2itt, tr2beta_s};

const ITERATIONS: usize = 100_000;

fn main() {
    println!("=== IAPWS-95 Performance Benchmark ===");
    println!("Iterations per test: {}\n", ITERATIONS);

    let t_c = 25.0;
    let rho = 997.0;

    let start = Instant::now();
    let mut sum = 0.0;
    for _ in 0..ITERATIONS { sum += tr2p(black_box(t_c), black_box(rho)); }
    let ms = start.elapsed().as_secs_f64() * 1000.0;
    println!("tr2p (pressure):          {:>10.2} ms ({:>8.0} ns/call)", ms, ms * 10000.0);

    sum = 0.0;
    let start = Instant::now();
    for _ in 0..ITERATIONS { sum += tr2u(black_box(t_c), black_box(rho)); }
    let ms = start.elapsed().as_secs_f64() * 1000.0;
    println!("tr2u (internal energy):   {:>10.2} ms ({:>8.0} ns/call)", ms, ms * 10000.0);

    sum = 0.0;
    let start = Instant::now();
    for _ in 0..ITERATIONS { sum += tr2h(black_box(t_c), black_box(rho)); }
    let ms = start.elapsed().as_secs_f64() * 1000.0;
    println!("tr2h (enthalpy):          {:>10.2} ms ({:>8.0} ns/call)", ms, ms * 10000.0);

    sum = 0.0;
    let start = Instant::now();
    for _ in 0..ITERATIONS { sum += tr2s(black_box(t_c), black_box(rho)); }
    let ms = start.elapsed().as_secs_f64() * 1000.0;
    println!("tr2s (entropy):           {:>10.2} ms ({:>8.0} ns/call)", ms, ms * 10000.0);

    sum = 0.0;
    let start = Instant::now();
    for _ in 0..ITERATIONS { sum += tr2cv(black_box(t_c), black_box(rho)); }
    let ms = start.elapsed().as_secs_f64() * 1000.0;
    println!("tr2cv (cv):               {:>10.2} ms ({:>8.0} ns/call)", ms, ms * 10000.0);

    sum = 0.0;
    let start = Instant::now();
    for _ in 0..ITERATIONS { sum += tr2cp(black_box(t_c), black_box(rho)); }
    let ms = start.elapsed().as_secs_f64() * 1000.0;
    println!("tr2cp (cp):               {:>10.2} ms ({:>8.0} ns/call)", ms, ms * 10000.0);

    sum = 0.0;
    let start = Instant::now();
    for _ in 0..ITERATIONS { sum += tr2w(black_box(t_c), black_box(rho)); }
    let ms = start.elapsed().as_secs_f64() * 1000.0;
    println!("tr2w (speed of sound):    {:>10.2} ms ({:>8.0} ns/call)", ms, ms * 10000.0);

    sum = 0.0;
    let start = Instant::now();
    for _ in 0..ITERATIONS { sum += tr2jt(black_box(t_c), black_box(rho)); }
    let ms = start.elapsed().as_secs_f64() * 1000.0;
    println!("tr2jt (Joule-Thomson):    {:>10.2} ms ({:>8.0} ns/call)", ms, ms * 10000.0);

    sum = 0.0;
    let start = Instant::now();
    for _ in 0..ITERATIONS { sum += tr2itt(black_box(t_c), black_box(rho)); }
    let ms = start.elapsed().as_secs_f64() * 1000.0;
    println!("tr2itt (isoth. throttle): {:>10.2} ms ({:>8.0} ns/call)", ms, ms * 10000.0);

    sum = 0.0;
    let start = Instant::now();
    for _ in 0..ITERATIONS { sum += tr2beta_s(black_box(t_c), black_box(rho)); }
    let ms = start.elapsed().as_secs_f64() * 1000.0;
    println!("tr2beta_s (isentropic):   {:>10.2} ms ({:>8.0} ns/call)", ms, ms * 10000.0);

    println!("\n=== Steam Properties (500C, 5 kg/m3) ===");
    let t_c = 500.0;
    let rho = 5.0;

    sum = 0.0;
    let start = Instant::now();
    for _ in 0..ITERATIONS {
        sum += tr2p(black_box(t_c), black_box(rho));
        sum += tr2h(black_box(t_c), black_box(rho));
        sum += tr2cp(black_box(t_c), black_box(rho));
        sum += tr2w(black_box(t_c), black_box(rho));
    }
    let ms = start.elapsed().as_secs_f64() * 1000.0;
    println!("Steam (p,h,cp,w) x 4:     {:>10.2} ms ({:>8.0} ns/call)", ms, ms * 10000.0);

    println!("\n=== All Properties x 5 States ===");
    let states = [
        (25.0, 997.0),
        (100.0, 0.6),
        (200.0, 10.0),
        (350.0, 950.0),
        (500.0, 5.0),
    ];

    sum = 0.0;
    let start = Instant::now();
    for _ in 0..ITERATIONS {
        for (t, r) in states {
            sum += tr2p(black_box(t), black_box(r));
            sum += tr2u(black_box(t), black_box(r));
            sum += tr2h(black_box(t), black_box(r));
            sum += tr2s(black_box(t), black_box(r));
            sum += tr2cv(black_box(t), black_box(r));
            sum += tr2cp(black_box(t), black_box(r));
            sum += tr2w(black_box(t), black_box(r));
            sum += tr2jt(black_box(t), black_box(r));
            sum += tr2itt(black_box(t), black_box(r));
            sum += tr2beta_s(black_box(t), black_box(r));
        }
    }
    let ms = start.elapsed().as_secs_f64() * 1000.0;
    println!("50 property calcs x 100K: {:>10.2} ms ({:>8.0} ns/state)", ms, ms * 10000.0);
}
