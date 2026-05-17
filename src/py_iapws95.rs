use pyo3::prelude::*;

use crate::iapws95::*;
use crate::iapws95_saturation::calc_saturation_properties;

#[pyfunction]
fn pt2h(p: f64, t_c: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    let rho = solve_density(p, t_k)?;
    Ok(calc_enthalpy(t_k, rho))
}

#[pyfunction]
fn pt2s(p: f64, t_c: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    let rho = solve_density(p, t_k)?;
    Ok(calc_entropy(t_k, rho))
}

#[pyfunction]
fn tx2h(t_c: f64, x: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    let sat = calc_saturation_properties(t_k)
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Temperature {} C is out of valid saturation range", t_c)
        ))?;
    Ok(sat.h_l + x * (sat.h_v - sat.h_l))
}

#[pyfunction]
fn tx2s(t_c: f64, x: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    let sat = calc_saturation_properties(t_k)
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Temperature {} C is out of valid saturation range", t_c)
        ))?;
    Ok(sat.s_l + x * (sat.s_v - sat.s_l))
}

#[pyfunction]
fn saturation_properties(t_c: f64) -> PyResult<(f64, f64, f64, f64, f64, f64, f64)> {
    let t_k = t_c + 273.15;
    let sat = calc_saturation_properties(t_k)
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Temperature {} C is out of valid saturation range", t_c)
        ))?;
    Ok((sat.p_sat, sat.rho_l, sat.rho_v, sat.h_l, sat.h_v, sat.s_l, sat.s_v))
}

fn solve_density(p: f64, t_k: f64) -> PyResult<f64> {
    if p <= 0.0 || t_k <= 0.0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid pressure or temperature"));
    }

    let rho_ideal = p * 1000.0 / (IAPWS95_R * t_k);

    let guesses = if rho_ideal < IAPWS95_RHOCRIT * 0.5 {
        vec![rho_ideal.max(1e-6), rho_ideal * 0.5, rho_ideal * 2.0, IAPWS95_RHOCRIT]
    } else {
        vec![IAPWS95_RHOCRIT, rho_ideal, rho_ideal * 0.5, rho_ideal * 2.0]
    };

    for rho_init in guesses {
        if let Ok(rho) = solve_density_newton(p, t_k, rho_init) {
            return Ok(rho);
        }
    }

    Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
        format!("Failed to solve density for p={} MPa, T={} K", p, t_k)
    ))
}

fn solve_density_newton(p: f64, t_k: f64, rho_init: f64) -> PyResult<f64> {
    let mut rho = rho_init.max(1e-8);

    for _ in 0..200 {
        let p_calc = calc_pressure(t_k, rho);
        let f = p_calc - p;
        if f.abs() < 1e-6 {
            return Ok(rho);
        }

        let drho = (rho * 1e-8).max(1e-10);
        let p_plus = calc_pressure(t_k, rho + drho);
        let df_drho = (p_plus - p_calc) / drho;

        if df_drho.abs() < 1e-20 {
            break;
        }

        let delta_rho = -f / df_drho;
        let damping = if delta_rho.abs() > rho * 0.5 { 0.1 } else { 0.5 };
        rho = (rho + damping * delta_rho).max(1e-8);
    }

    let p_calc = calc_pressure(t_k, rho);
    if (p_calc - p).abs() < 1e-3 {
        Ok(rho)
    } else {
        Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Newton iteration failed"))
    }
}

#[pymodule]
fn iapws95(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pt2h, m)?)?;
    m.add_function(wrap_pyfunction!(pt2s, m)?)?;
    m.add_function(wrap_pyfunction!(tx2h, m)?)?;
    m.add_function(wrap_pyfunction!(tx2s, m)?)?;
    m.add_function(wrap_pyfunction!(saturation_properties, m)?)?;
    Ok(())
}
