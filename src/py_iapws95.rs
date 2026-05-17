use pyo3::prelude::*;

use crate::iapws95::*;
use crate::iapws95_saturation::calc_saturation_properties;

// ==========================================================================
// Functions for (p,T) → property calculations using numerical inversion
// ==========================================================================

#[pyfunction]
fn pt2h(p: f64, t_c: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    let rho = solve_density(p, t_k).ok_or_else(|| 
        PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Failed to solve density for p={} MPa, T={} K", p, t_k)
        )
    )?;
    Ok(calc_enthalpy(t_k, rho))
}

#[pyfunction]
fn pt2s(p: f64, t_c: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    let rho = solve_density(p, t_k).ok_or_else(|| 
        PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Failed to solve density for p={} MPa, T={} K", p, t_k)
        )
    )?;
    Ok(calc_entropy(t_k, rho))
}

// ==========================================================================
// Functions for (T,rho) → property calculations - direct computation
// ==========================================================================

#[pyfunction]
fn tr2p(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_pressure(t_k, rho))
}

#[pyfunction]
fn tr2u(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_internal_energy(t_k, rho))
}

#[pyfunction]
fn tr2h(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_enthalpy(t_k, rho))
}

#[pyfunction]
fn tr2s(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_entropy(t_k, rho))
}

#[pyfunction]
fn tr2cv(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_cv(t_k, rho))
}

#[pyfunction]
fn tr2cp(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_cp(t_k, rho))
}

#[pyfunction]
fn tr2w(t_c: f64, rho: f64) -> PyResult<f64> {
    let t_k = t_c + 273.15;
    Ok(calc_speed_of_sound(t_k, rho))
}

// ==========================================================================
// Functions for (T,x) → property calculations in two-phase region
// ==========================================================================

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

// ==========================================================================
// Saturation properties
// ==========================================================================

#[pyfunction]
fn saturation_properties(t_c: f64) -> PyResult<(f64, f64, f64, f64, f64, f64, f64)> {
    let t_k = t_c + 273.15;
    let sat = calc_saturation_properties(t_k)
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Temperature {} C is out of valid saturation range", t_c)
        ))?;
    Ok((sat.p_sat, sat.rho_l, sat.rho_v, sat.h_l, sat.h_v, sat.s_l, sat.s_v))
}

#[pymodule]
fn iapws95(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // (p,T) → property calculations
    m.add_function(wrap_pyfunction!(pt2h, m)?)?;
    m.add_function(wrap_pyfunction!(pt2s, m)?)?;
    
    // (T,rho) → property calculations
    m.add_function(wrap_pyfunction!(tr2p, m)?)?;
    m.add_function(wrap_pyfunction!(tr2u, m)?)?;
    m.add_function(wrap_pyfunction!(tr2h, m)?)?;
    m.add_function(wrap_pyfunction!(tr2s, m)?)?;
    m.add_function(wrap_pyfunction!(tr2cv, m)?)?;
    m.add_function(wrap_pyfunction!(tr2cp, m)?)?;
    m.add_function(wrap_pyfunction!(tr2w, m)?)?;
    
    // (T,x) → property calculations
    m.add_function(wrap_pyfunction!(tx2h, m)?)?;
    m.add_function(wrap_pyfunction!(tx2s, m)?)?;
    
    // Saturation properties
    m.add_function(wrap_pyfunction!(saturation_properties, m)?)?;
    
    Ok(())
}
