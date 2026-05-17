use iapws95::iapws95::{calc_entropy, calc_enthalpy, IAPWS95_TCRIT, solve_density};
use iapws95::iapws95_saturation::calc_saturation_properties;

use plotters::prelude::*;

const P_MIN: f64 = 611.657e-6;
const S_MIN: f64 = 0.0;
const S_MAX: f64 = 12.5;
const H_MIN: f64 = 0.0;
const H_MAX: f64 = 4300.0;

fn generate_isotherm_data(t_c: f64, pressures: &[f64]) -> Vec<(f64, f64)> {
    let t_k = t_c + 273.15;
    let mut points = Vec::new();

    for &p in pressures {
        if let Some(rho) = solve_density(p, t_k) {
            let h = calc_enthalpy(t_k, rho);
            let s = calc_entropy(t_k, rho);
            if s >= S_MIN && s <= S_MAX && h >= H_MIN && h <= H_MAX {
                points.push((s, h));
            }
        }
    }

    points
}

fn generate_isobar_data(p: f64, temperatures_c: &[f64]) -> Vec<(f64, f64)> {
    let mut points = Vec::new();

    for &t_c in temperatures_c {
        let t_k = t_c + 273.15;
        if let Some(rho) = solve_density(p, t_k) {
            let h = calc_enthalpy(t_k, rho);
            let s = calc_entropy(t_k, rho);
            if s >= S_MIN && s <= S_MAX && h >= H_MIN && h <= H_MAX {
                points.push((s, h));
            }
        }
    }

    points
}

fn generate_saturation_data() -> (Vec<(f64, f64)>, Vec<(f64, f64)>) {
    let mut liquid_points = Vec::new();
    let mut vapor_points = Vec::new();

    let t_min = 273.16;
    let t_max = IAPWS95_TCRIT;
    let num_points = 200;

    for i in 0..=num_points {
        let t = t_min + (t_max - t_min) * i as f64 / num_points as f64;
        if let Some(sat) = calc_saturation_properties(t) {
            if sat.s_l >= S_MIN && sat.s_l <= S_MAX && sat.h_l >= H_MIN && sat.h_l <= H_MAX {
                liquid_points.push((sat.s_l, sat.h_l));
            }
            if sat.s_v >= S_MIN && sat.s_v <= S_MAX && sat.h_v >= H_MIN && sat.h_v <= H_MAX {
                vapor_points.push((sat.s_v, sat.h_v));
            }
        }
    }

    (liquid_points, vapor_points)
}

fn generate_isoquality_data(x: f64) -> Vec<(f64, f64)> {
    let mut points = Vec::new();

    let t_min = 273.16;
    let t_max = IAPWS95_TCRIT - 0.1;
    let num_points = 100;

    for i in 0..=num_points {
        let t = t_min + (t_max - t_min) * i as f64 / num_points as f64;
        if let Some(sat) = calc_saturation_properties(t) {
            let h = sat.h_l + x * (sat.h_v - sat.h_l);
            let s = sat.s_l + x * (sat.s_v - sat.s_l);
            if s >= S_MIN && s <= S_MAX && h >= H_MIN && h <= H_MAX {
                points.push((s, h));
            }
        }
    }

    points
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== IAPWS-95 Mollier (H-S) Diagram Generator ===\n");

    let root = BitMapBackend::new("mollier_diagram.png", (1200, 900)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("H-S (Mollier) Diagram", ("sans-serif", 30).into_font())
        .margin(20)
        .x_label_area_size(60)
        .y_label_area_size(80)
        .build_cartesian_2d(S_MIN..S_MAX, H_MIN..H_MAX)?;

    chart
        .configure_mesh()
        .x_desc("s, kJ/(kg·K)")
        .y_desc("h, kJ/kg")
        .x_labels(13)
        .y_labels(10)
        .label_style(("sans-serif", 12))
        .draw()?;

    let isotherms_c = vec![0.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0];
    let isobars_mpa = vec![
        P_MIN,
        0.001,
        0.01,
        0.1,
        1.0,
        10.0,
        20.0,
        50.0,
        100.0,
    ];

    println!("Calculating isotherm lines...");
    for &t in &isotherms_c {
        let points = generate_isotherm_data(t, &isobars_mpa);
        if points.len() > 1 {
            chart
                .draw_series(LineSeries::new(
                    points.into_iter().map(|(s, h)| (s, h)),
                    &GREEN,
                ))?
                .label(format!("{}°C", t))
                .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));
        }
    }

    println!("Calculating isobar lines...");
    let isobar_temps_c: Vec<f64> = (0..=800).step_by(50).map(|t| t as f64).collect();
    for &p in &isobars_mpa {
        let points = generate_isobar_data(p, &isobar_temps_c);
        if points.len() > 1 {
            chart
                .draw_series(LineSeries::new(
                    points.into_iter().map(|(s, h)| (s, h)),
                    &BLUE,
                ))?
                .label(format!("{:.4} MPa", p))
                .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
        }
    }

    println!("Calculating saturation lines...");
    let (liquid_sat, vapor_sat) = generate_saturation_data();

    if !liquid_sat.is_empty() {
        chart
            .draw_series(LineSeries::new(
                liquid_sat.into_iter().map(|(s, h)| (s, h)),
                ShapeStyle::from(&RED).stroke_width(2),
            ))?
            .label("x=0 (saturated liquid)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    }

    if !vapor_sat.is_empty() {
        chart
            .draw_series(LineSeries::new(
                vapor_sat.into_iter().map(|(s, h)| (s, h)),
                ShapeStyle::from(&RED).stroke_width(2),
            ))?
            .label("x=1 (saturated vapor)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    }

    println!("Calculating isoquality lines...");
    let isoqualities: Vec<f64> = (1..=9).map(|i| i as f64 / 10.0).collect();
    for x in &isoqualities {
        let points = generate_isoquality_data(*x);
        if points.len() > 1 {
            chart.draw_series(LineSeries::new(
                points.into_iter().map(|(s, h)| (s, h)),
                ShapeStyle::from(&RGBColor(255, 100, 100)).stroke_width(1),
            ))?;
        }
    }

    chart
        .configure_series_labels()
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    println!("\nMollier diagram saved to: mollier_diagram.png");
    println!("Diagram range: s = {} to {} kJ/(kg·K), h = {} to {} kJ/kg", S_MIN, S_MAX, H_MIN, H_MAX);

    Ok(())
}
